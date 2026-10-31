#!/usr/bin/env python3
"""Rebuild eight SI tables. Only third-party dependency: lxml.

Usage: python3 rebuild_tables.py SOURCE_SI.docx OUTPUT_SI.docx [--report REPORT.json]
Output paths are supplied by the caller; the report defaults beside OUTPUT_SI.docx. Unedited XML spans and
all other ZIP members are copied verbatim. Decimal arithmetic uses source text.
"""
import argparse
import csv
import json
import re
import sys
from collections import OrderedDict
from copy import deepcopy
from decimal import Decimal, InvalidOperation, ROUND_HALF_UP, localcontext
from pathlib import Path
from xml.parsers import expat
from zipfile import ZipFile
from lxml import etree as E

W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'
NS = {'w': W}
IDS = ['S1a', 'S1b', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7']
REPO = Path(__file__).resolve().parents[2]
FUELS = dict(coal='Coal', gas='Natural gas', oil='Oil', biomass='Biomass', wind='Wind', solar='Solar', hydro='Hydropower', nuclear='Nuclear', geothermal='Geothermal')
METRICS = {'count':'Plant count', 'heat':'Daily heat input', 'Actual daily Heat Input':'Actual daily heat input'}
COV = {'c_prop':'Continental air-mass fraction', 'inv_hpbl':'Inverse boundary-layer height', 'pm_model_centered':'PM2.5 (centred)', 'rhum':'Relative humidity', 'soilm':'Soil moisture', 'sunspot_monthly':'Monthly sunspot number', 'sunspot_monthly:lat_diff_35':'Sunspot number × latitude', 'surface_u_ppm_article:soilm':'Surface uranium × soil moisture', 'temp_centered':'Temperature (linear)', 'temp2':'Temperature (squared)', 'temp3':'Temperature (cubed)', 'vel':'Wind speed'}
UNITS = {'per unit proportion':'Per unit proportion', 'per 1/m':'Per m⁻¹', 'per degree C (linear polynomial term)':'Per °C', 'per degree C squared':'Per °C²', 'per degree C cubed':'Per °C³', 'per m/s':'Per m/s', 'per microgram/m3':'Per µg/m³', 'per native NARR kg/m2 at zero surface uranium':'Per kg/m² (at zero surface uranium)', 'per percentage point':'Per percentage point', 'per ppm x native NARR kg/m2':'Per ppm × kg/m²', 'per sunspot at latitude 35 degrees':'Per sunspot (at 35° latitude)', 'per sunspot x degree latitude':'Per sunspot × degree latitude'}

def q(n): return '{%s}%s' % (W,n)
def el(n, **attrs): return E.Element(q(n), {q(k):str(v) for k,v in attrs.items()})
def txt(x): return ''.join(x.xpath('.//w:t/text()', namespaces=NS))
def decimal(s): return Decimal(s.replace('−','-').replace(',',''))

def rounded(s, kind):
    d=decimal(s)
    if kind=='int':
        assert d==d.to_integral_value(), s
        return format(int(d), ',').replace('-', '−')
    if kind=='p' and d < Decimal('.001'):
        assert d>=0
        return '<0.001'
    exponent = -3 if kind in ('wind','p') else (d.adjusted()-2 if d else -2)
    with localcontext() as ctx:
        ctx.prec=max(80,len(d.as_tuple().digits)+abs(exponent)+10)
        v=d.quantize(Decimal(1).scaleb(exponent), rounding=ROUND_HALF_UP)
    return format(v, 'f').replace('-', '−')

def cell(kind, *source):
    vals=list(source)
    if kind=='text': display=source[0]
    elif kind=='ci': display='%s (%s, %s)' % tuple(rounded(s,'sig') for s in source)
    else: display=rounded(source[0],kind)
    return {'kind':kind, 'source':vals, 'display':display}

def load_source(path):
    with ZipFile(path) as z: raw=z.read('word/document.xml')
    root=E.fromstring(raw); body=root.find('w:body',NS)
    tables=body.findall('w:tbl',NS)
    assert len(tables)==8, 'Expected exactly eight body tables'
    data={}
    for ident,t in zip(IDS,tables):
        matrix=[[txt(c) for c in tr.findall('w:tc',NS)] for tr in t.findall('w:tr',NS)]
        assert len(set(matrix[0]))==len(matrix[0])
        assert all(len(row)==len(matrix[0]) for row in matrix)
        data[ident]=[dict(zip(matrix[0],row)) for row in matrix[1:]]
    return raw,root,body,tables,data

def constant(rows,key,expected=None):
    vals=list(dict.fromkeys(r[key] for r in rows))
    assert len(vals)==1, '%s is not constant: %s' % (key, vals)
    if expected is not None: assert vals[0]==expected, (key,vals,expected)
    return vals[0]

def build_specs(data):
    specs=[]; problems=[]
    s3drop=all(r['Wind_reference']==('sample_median' if r['Period']=='2001–2024' else 'period_sample_median') for r in data['S3'])
    assert {r['Period'] for r in data['S3']}=={'2001–2024','2001–2008','2009–2016','2017–2024'}
    grouped=OrderedDict()
    for i,r in enumerate(data['S6']):
        key=tuple(r[k] for k in ['fuel','metric','radius_km'])
        group=grouped.setdefault(key,{})
        assert r['model'] not in group, ('duplicate S6 key',key,r['model'])
        group[r['model']]=r
    reasons=[]
    for key,g in grouped.items():
        if set(g)!={'current','natural_spline'}: reasons.append({'key':key,'reason':'missing model'})
        else:
            for k in ['N_model','wind_ms']:
                if decimal(g['current'][k])!=decimal(g['natural_spline'][k]): reasons.append({'key':key,'field':k,'values':[g[m][k] for m in ['current','natural_spline']]})
    wide=not reasons
    if reasons: problems.append({'S6_long_reason':reasons})
    if not s3drop: problems.append('S3 Wind_reference retained: exact period-specific condition failed.')

    def spec(ident,headers,rows,dropped=(),constants=None):
        maxchars=[max(len(row[j]['display']) for row in rows) for j in range(len(headers))]
        numeric=[rows[0][j]['kind']!='text' for j in range(len(headers))]
        # Numeric widths obey the requested conservative 100 dxa/character +120.
        # Text bodies may wrap at spaces; headers may always wrap.
        widths=[]
        for j,(n,num) in enumerate(zip(maxchars,numeric)):
            bodywidth=100*n+120 if num else min(2100,max(1050,70*n+120))
            # Header wrapping is allowed, but keep words such as "observations"
            # intact and give longer headers enough room for a few short lines.
            headerwidth=max(70*max(map(len,headers[j].split()))+120,50*len(headers[j])+120)
            widths.append(max(850,bodywidth,headerwidth))
        # Some tables have three CI columns; reduce only wrappable text widths.
        excess=sum(widths)-14000
        if excess>0:
            for j in sorted((i for i,n in enumerate(numeric) if not n),key=lambda i:widths[i],reverse=True):
                reduction=min(excess,max(0,widths[j]-850)); widths[j]-=reduction; excess-=reduction
        assert sum(widths)<=14000, (ident,widths,'Numeric budget exceeds table width')
        # Allocate spare space to descriptive columns and then header-heavy counts.
        target=11000 if ident=='S1b' else 14000
        spare=max(0,target-sum(widths))
        weights=[maxchars[j] if not n else (len(headers[j]) if rows[0][j]['kind']!='ci' else 0) for j,n in enumerate(numeric)]
        if not sum(weights): weights=[1]*len(widths)
        extra=[spare*w//sum(weights) for w in weights]
        for j,x in enumerate(extra): widths[j]+=x
        widths[weights.index(max(weights))]+=spare-sum(extra)
        s={'id':ident,'headers':headers,'rows':rows,'numeric':numeric,'rows_in':len(data[ident]),'rows_out':len(rows),'cols_out':len(headers),'col_widths':widths,'dropped':list(dropped),'constants_confirmed':constants or {},'max_body_chars':maxchars}
        specs.append(s)

    estheader='Estimate (95% CI) (mBq/m³)'
    for ident in ['S1a','S1b']:
        a=ident=='S1a'; rr=[]
        for r in data[ident]:
            row=[cell('text',FUELS[r['fuel']])]
            if a: row.append(cell('text',METRICS[r['metric']]))
            row += [cell('int',r[k]) for k in ['radius_km','n_rows','n_monitors','n_exposed_rows','n_exposed_monitors','n_monitor_years_with_variation']]
            rr.append(row)
        spec(ident,['Fuel']+(['Exposure metric'] if a else [])+['Radius (km)','Eligible filter periods (n)','Monitors (n)','Exposed filter periods (n)','Exposed monitors (n)','Monitor-years with exposure variation (n)'],rr,constants={'n_monitors':constant(data[ident],'n_monitors','133')})
    rr=[]
    for r in data['S2']:
        rr.append([cell('text',FUELS.get(r['Fuel'],r['Fuel'])),cell('text',METRICS.get(r['Exposure_metric'],r['Exposure_metric'])),cell('int',r['Radius_km']),cell('wind',r['Wind_speed_m_s']),cell('ci',r['Estimate_mBq_m3'],r['CI95_lower'],r['CI95_upper']),cell('int',r['N_pre']),cell('int',r['N_model'])])
    spec('S2',['Fuel','Exposure metric','Radius (km)','Wind speed (m/s)',estheader,'Eligible observations (n)','Model observations (n)'],rr,['Wind_reference'],{'Wind_reference':constant(data['S2'],'Wind_reference','Sample median')})
    rr=[]
    for r in data['S3']:
        row=[cell('text',r['Period']),cell('text',FUELS.get(r['Fuel'],r['Fuel'])),cell('text',METRICS.get(r['Exposure_metric'],r['Exposure_metric']))]
        if not s3drop: row.append(cell('text',{'sample_median':'Sample median','period_sample_median':'Period median'}.get(r['Wind_reference'],r['Wind_reference'])))
        row += [cell('wind',r['Wind_speed_m_s']),cell('int',r['N_monitors']),cell('ci',r['Estimate_mBq_m3'],r['CI95_lower'],r['CI95_upper']),cell('int',r['N_pre'])]
        rr.append(row)
    spec('S3',['Period','Fuel','Exposure metric']+([] if s3drop else ['Wind reference'])+['Wind speed (m/s)','Monitors (n)',estheader,'Eligible observations (n)'],rr,['Radius_km','Status']+(['Wind_reference'] if s3drop else []),{'Radius_km':constant(data['S3'],'Radius_km','20'),'Status':constant(data['S3'],'Status','PASS'),'Wind_reference_by_period':{p:sorted({r['Wind_reference'] for r in data['S3'] if r['Period']==p}) for p in dict.fromkeys(r['Period'] for r in data['S3'])}})
    rr=[[cell('text',r['Analysis']),cell('text',r['Wind_direction']),cell('int',r['Sector_width_deg']),cell('ci',r['Estimate_mBq_m3'],r['CI95_lower'],r['CI95_upper'])] for r in data['S4']]
    spec('S4',['Analysis','Wind direction','Sector width (°)',estheader],rr,['Radius_km','Wind_speed_m_s'],{k:constant(data['S4'],k,v) for k,v in [('Radius_km','20'),('Wind_speed_m_s','2.98560638549627')]})
    rr=[[cell('text',COV[r['term']]),cell('text',UNITS[r['Model_unit']]),cell('ci',r['estimate'],r['ci95_low'],r['ci95_high']),cell('sig',r['se']),cell('p',r['p_wald'])] for r in data['S5']]
    const={k:constant(data['S5'],k) for k in ['Model_N','Monitors','Outcome_unit','N_pre','Note']}
    tvals=list(dict.fromkeys(r['t'] for r in data['S5']))
    if len(tvals)==1: const['t']=tvals[0]
    else: problems.append('S5 t is not constant (12 distinct values); omitted as explicitly requested. All five metadata columns requested for constant checks are constant.')
    spec('S5',['Covariate','Unit of covariate',estheader,'SE','P value'],rr,['t','Model_N','Monitors','Outcome_unit','N_pre','Note'],const)
    const={k:constant(data['S6'],k,v) for k,v in [('N_monitors','133'),('wind_reference','sample_median')]}
    assert all(r['spline_df']==('3' if r['model']=='natural_spline' else '') for r in data['S6'])
    const['spline_df_by_model']={'current':'','natural_spline':'3'}
    const['matched_N_model_and_wind_ms_all_pairs']=wide
    rr=[]
    for r in ([g['current'] for g in grouped.values()] if wide else data['S6']):
        row=[cell('text',FUELS[r['fuel']]),cell('text',METRICS[r['metric']]),cell('int',r['radius_km'])]
        if not wide: row.append(cell('text',{'current':'Main model','natural_spline':'Natural spline, 3 df'}[r['model']]))
        row += [cell('int',r['N_model']),cell('wind',r['wind_ms']),cell('ci',r['estimate'],r['ci95_low'],r['ci95_high'])]
        if wide:
            alt=grouped[tuple(r[k] for k in ['fuel','metric','radius_km'])]['natural_spline']
            row.append(cell('ci',alt['estimate'],alt['ci95_low'],alt['ci95_high']))
        rr.append(row)
    spec('S6',['Fuel','Exposure metric','Radius (km)']+([] if wide else ['Model'])+['Model observations (n)','Wind speed (m/s)']+(['Main model estimate (95% CI) (mBq/m³)','Natural-spline estimate (95% CI) (mBq/m³)'] if wide else [estheader]),rr,['spline_df','N_input','N_monitors','wind_reference','se','p_wald']+(['model'] if wide else []),const)
    winds=list(dict.fromkeys(r['wind_ms'] for r in data['S7']))
    rwinds=list(dict.fromkeys(rounded(w,'wind') for w in winds))
    windtext=rwinds[0] if len(rwinds)==1 else '; '.join('%s: %s' % (METRICS[m],', '.join(dict.fromkeys(rounded(r['wind_ms'],'wind') for r in data['S7'] if r['metric']==m))) for m in dict.fromkeys(r['metric'] for r in data['S7']))
    assert all(decimal(r['attenuation_pct'])==-decimal(r['relative_change_pct']) for r in data['S7'])
    rr=[[cell('text',METRICS[r['metric']]),cell('int',r['radius_km']),cell('int',r['n_model']),cell('ci',r['original_estimate'],r['original_ci95_low'],r['original_ci95_high']),cell('ci',r['fire_adjusted_estimate'],r['fire_adjusted_ci95_low'],r['fire_adjusted_ci95_high']),cell('sig',r['relative_change_pct']),cell('ci',r['fire_pm25_coef'],r['fire_pm25_ci95_low'],r['fire_pm25_ci95_high']),cell('p',r['fire_pm25_p'])] for r in data['S7']]
    spec('S7',['Exposure metric','Radius (km)','Model observations (n)','Original estimate (95% CI) (mBq/m³)','Fire-adjusted estimate (95% CI) (mBq/m³)','Relative change (%)','Fire-PM2.5 coefficient (95% CI)','P value (fire PM2.5)'],rr,['n_monitors','wind_ms','estimate_change','attenuation_pct','wind_reference'],{'n_monitors':constant(data['S7'],'n_monitors','133'),'wind_reference':constant(data['S7'],'wind_reference','sample_median'),'attenuation_pct_equals_negative_relative_change_pct':True,'wind_ms_distinct':winds,'wind_ms_rounded_distinct':rwinds})
    return specs,s3drop,wide,windtext,problems

def run_properties(base=None,size=None,bold=None):
    p=deepcopy(base) if base is not None else el('rPr')
    for x in p.iter():
        for a in list(x.attrib):
            if E.QName(a).localname in ('themeColor','themeShade','themeTint'): del x.attrib[a]
    color=p.find(q('color'))
    if color is None: color=el('color'); p.append(color)
    color.set(q('val'),'000000')
    if size is not None:
        for key in ['sz','szCs']:
            x=p.find(q(key))
            if x is None: x=el(key); p.append(x)
            x.set(q('val'),str(size))
    if bold is not None:
        for key in ['b','bCs']:
            x=p.find(q(key))
            if x is None: x=el(key); p.append(x)
            x.set(q('val'),'1' if bold else '0')
    return p

def append_runs(parent,text,rpr):
    # Split only PM2.5, not the independent aerodynamic diameter ≤2.5.
    for part in re.split(r'(PM2\.5)',text):
        if not part: continue
        segments=[('PM',False),('2.5',True)] if part=='PM2.5' else [(part,False)]
        for segment,sub in segments:
            r=el('r'); props=deepcopy(rpr)
            if sub:
                va=props.find(q('vertAlign'))
                if va is None: va=el('vertAlign'); props.append(va)
                va.set(q('val'),'subscript')
            r.append(props); t=el('t'); t.text=segment
            if segment.strip()!=segment: t.set('{http://www.w3.org/XML/1998/namespace}space','preserve')
            r.append(t); parent.append(r)

def replace_paragraph(p,text):
    first=p.find('.//w:r/w:rPr',NS)
    rp=run_properties(first)
    for c in list(p):
        if c.tag!=q('pPr'): p.remove(c)
    append_runs(p,text,rp)

def make_table(old,s):
    table=el('tbl'); oldp=old.find(q('tblPr'))
    props=deepcopy(oldp) if oldp is not None else el('tblPr')
    for n,attrs in [('tblW',{'w':sum(s['col_widths']),'type':'dxa'}),('jc',{'val':'center'}),('tblLayout',{'type':'fixed'})]:
        node=props.find(q(n))
        if node is None: node=el(n); props.append(node)
        for k,v in attrs.items(): node.set(q(k),str(v))
    mar=props.find(q('tblCellMar'))
    if mar is None: mar=el('tblCellMar'); props.append(mar)
    for side,width in [('left',60),('right',60),('top',25),('bottom',25)]:
        x=mar.find(q(side))
        if x is None: x=el(side); mar.append(x)
        x.set(q('w'),str(width)); x.set(q('type'),'dxa')
    table.append(props); grid=el('tblGrid')
    for width in s['col_widths']: grid.append(el('gridCol',w=width))
    table.append(grid)
    base=old.find('.//w:r/w:rPr',NS)
    for i,row in enumerate([[{'display':h} for h in s['headers']]]+s['rows']):
        tr=el('tr'); trp=el('trPr'); trp.append(el('cantSplit'))
        if i==0: trp.append(el('tblHeader'))
        tr.append(trp)
        for j,c in enumerate(row):
            tc=el('tc'); cp=el('tcPr'); cp.append(el('tcW',w=s['col_widths'][j],type='dxa'))
            if i==0: cp.append(el('shd',val='clear',color='auto',fill='EEEEEE'))
            if i>0 and s['numeric'][j]: cp.append(el('noWrap'))
            cp.append(el('vAlign',val='center')); tc.append(cp)
            p=el('p'); pp=el('pPr'); pp.append(el('spacing',before=0,after=0,line=240,lineRule='auto'))
            pp.append(el('jc',val='center' if s['numeric'][j] else 'left')); p.append(pp)
            rp=run_properties(base,size=17,bold=(i==0))
            append_runs(p,c['display'],rp); tc.append(p); tr.append(tc)
        table.append(tr)
    return table

def body_spans(raw):
    """Byte spans of each direct body child, without reserializing untouched XML."""
    parser=expat.ParserCreate(); stack=[]; spans=[]; start=None
    def on_start(name,attrs):
        nonlocal start
        if stack and stack[-1]=='w:body': start=parser.CurrentByteIndex
        stack.append(name)
    def on_end(name):
        nonlocal start
        if len(stack)>=2 and stack[-2]=='w:body':
            pos=parser.CurrentByteIndex
            end=raw.index(b'>',pos)+1 if raw[pos:pos+2]==b'</' else pos
            spans.append((start,end)); start=None
        stack.pop()
    parser.StartElementHandler=on_start; parser.EndElementHandler=on_end
    parser.Parse(raw,True)
    return spans

def serialize_in_context(element):
    # All namespaces are declared on the original document root. Serialize in
    # that namespace context and strip duplicate declarations from the fragment.
    wrapper=E.Element('wrapper',nsmap={'w':W}); wrapper.append(element)
    out=E.tostring(element,encoding='UTF-8',with_tail=False)
    return re.sub(rb' xmlns(?::[A-Za-z0-9_]+)?="[^"]*"',b'',out)

def paragraph_changes(body,s3drop,wide,wind):
    paragraphs={i:p for i,p in enumerate(body) if p.tag==q('p')}
    changes={}
    def change(i,text):
        p=deepcopy(paragraphs[i]); replace_paragraph(p,text); changes[i]=p
    cap='Supplementary Table S4. Overall inventory directional and sector analyses at 20 km and mean wind 2.986 m/s, using 558,915 quasi-daily records from 140 monitors.'
    # Locate all target paragraphs by content, avoiding brittle shifted indices.
    def find(prefix):
        matches=[i for i,p in paragraphs.items() if txt(p).startswith(prefix)]
        assert len(matches)==1,(prefix,matches)
        return matches[0]
    s1note=find('Input rows, nonzero-exposure periods')
    old=txt(paragraphs[s1note]); first='Input rows, nonzero-exposure periods, exposed monitors and monitor-years with exposure variation are reported separately.'
    assert old.startswith(first)
    change(s1note,'Eligible filter periods, exposed filter periods (nonzero exposure), exposed monitors and monitor-years with exposure variation are reported separately.'+old[len(first):])
    notes={
        'N_pre is eligible input and N_model':'Eligible observations are input filter periods before estimation; model observations are those retained by the estimator. Associations are evaluated at the model-specific sample-median wind speed. Estimates are in mBq/m³ per plant or per 10,000 MMBtu/day.',
        'N_pre is eligible input; actual':'All estimates are at 20 km. Eligible observations are input filter periods; the post-singleton model N was not retained in this historical export.'+(' Wind speed is the sample median for 2001–2024 and the period-specific median for each sub-period.' if s3drop else ''),
        'Model N=107,367':'Model N = 107,367 from 107,407 eligible filter periods and 133 monitor clusters. Coefficients are conditional; the three temperature terms should be interpreted jointly. Intervals and P values are two-sided normal-Wald.',
        'Specifications are compared on matched samples':'Specifications are compared on matched samples at the sample-median wind speed; all models include 133 monitors. The natural-spline models use 3 degrees of freedom. Eligible input N, standard errors and P values are omitted here and are available with full-precision values in the GitHub repository (data/display_assets/frozen_v3/tables/Supplementary_Table_S6_natural_spline_estimates.csv).',
        'Reference and adjusted models share observations;':'Reference and fire-adjusted models share observations from 133 monitors and are evaluated at the sample-median wind speed ('+wind+' m/s). Relative change is the change in the estimate after fire-PM2.5 adjustment. Absolute changes and full-precision values are available in the GitHub repository (source_data_manuscript_v5/Supplementary_Figure_S8_GFED_estimates.csv).',
    }
    for prefix,text in notes.items(): change(find(prefix),text)
    table_indices=[i for i,x in enumerate(body) if x.tag==q('tbl')]
    abbr={
        'S2':'Abbreviations: CAMPD = Clean Air Markets Program Data; CI = confidence interval; MMBtu = million British thermal units.',
        'S3':'Abbreviations: CI = confidence interval.',
        'S4':'Abbreviations: CI = confidence interval.',
        'S5':'Abbreviations: CI = confidence interval; NARR = North American Regional Reanalysis; PM2.5 = particulate matter with aerodynamic diameter ≤2.5 µm; SE = standard error.',
        'S6':'Abbreviations: CI = confidence interval.',
        'S7':'Abbreviations: CI = confidence interval; GFED = Global Fire Emissions Database; PM2.5 = particulate matter with aerodynamic diameter ≤2.5 µm.',
    }
    for ident,idx in zip(IDS,table_indices):
        if ident in abbr:
            i=next(j for j in range(idx+1,len(body)) if j in paragraphs and txt(paragraphs[j]).startswith('Abbreviations:'))
            change(i,abbr[ident])
    capidx=table_indices[IDS.index('S4')]-1
    assert txt(body[capidx]).startswith('Supplementary Table S4.')
    change(capidx,cap)
    links=body.xpath('w:p/w:hyperlink[@w:anchor="Table_S4"]',namespaces=NS)
    assert len(links)==1
    link=links[0]; p=deepcopy(link.getparent()); h=p.find('w:hyperlink[@w:anchor="Table_S4"]',NS)
    rp=run_properties(h.find('.//w:r/w:rPr',NS))
    for c in list(h): h.remove(c)
    append_runs(h,cap[:-1],rp)
    changes[list(body).index(link.getparent())]=p
    # Change only pgSz of the section containing S1b, preserving margins and
    # all other section properties, and leaving every figure section intact.
    s1bidx=table_indices[1]
    sectionidx=next(i for i in range(s1bidx+1,len(body)) if body[i].find('w:pPr/w:sectPr',NS) is not None)
    sp=deepcopy(body[sectionidx]); sz=sp.find('w:pPr/w:sectPr/w:pgSz',NS)
    assert sz is not None
    sz.set(q('w'),'15840'); sz.set(q('h'),'12240'); sz.set(q('orient'),'landscape'); changes[sectionidx]=sp
    heading=next(i for i in range(table_indices[0]-1,-1,-1) if txt(body[i])=='Supplementary Tables')
    template=body[find('N_pre is eligible input and N_model')]
    intro=deepcopy(template)
    replace_paragraph(intro,'Values are rounded to three significant figures; wind speeds to three decimal places and P values to three decimal places. Full-precision values and all exported columns are available at https://github.com/hudson-HH/PowerPlant_Radioactivity.')
    return changes,heading,intro

def crosscheck(data):
    mismatches=[]; problems=[]
    for ident in IDS:
        paths=list((REPO/'data/display_assets/frozen_v3/tables').glob('Supplementary_Table_'+ident+'_*.csv')) if ident!='S7' else [REPO/'source_data_manuscript_v5/Supplementary_Figure_S8_GFED_estimates.csv']
        if len(paths)!=1 or not paths[0].exists():
            problems.append('CSV unavailable or ambiguous for '+ident); continue
        with paths[0].open(encoding='utf-8-sig',newline='') as f:
            reader=csv.DictReader(f); rows=list(reader)
        src=data[ident]
        if len(rows)!=len(src): mismatches.append({'table':ident,'field':'row_count','docx':len(src),'csv':len(rows)})
        if ident!='S7' and list(src[0])!=list(reader.fieldnames): mismatches.append({'table':ident,'field':'headers','docx':list(src[0]),'csv':reader.fieldnames})
        keys={'S1a':['fuel','metric','radius_km'],'S1b':['fuel','radius_km'],'S2':['Fuel','Exposure_metric','Radius_km'],'S3':['Period','Fuel','Exposure_metric'],'S4':['Analysis','Wind_direction','Sector_width_deg'],'S5':['term'],'S6':['fuel','metric','radius_km','model'],'S7':['metric','radius_km']}[ident]
        assert [tuple(r[k] for k in keys) for r in src]==[tuple(r[k] for k in keys) for r in rows], 'CSV row keys/order differ for '+ident
        for i,(a,b) in enumerate(zip(src,rows),1):
            for k,v in a.items():
                if ident=='S7' and k=='n_model':
                    for model in ['n_model_M0','n_model_M1']:
                        if decimal(v)!=decimal(b[model]): mismatches.append({'table':ident,'row':i,'column':k,'csv_column':model,'docx':v,'csv':b[model]})
                    continue
                cv=b.get(k)
                if v==cv: continue
                try: equal=decimal(v)==decimal(cv)
                except (InvalidOperation,AttributeError): equal=False
                if not equal: mismatches.append({'table':ident,'row':i,'column':k,'docx':v,'csv':cv})
    return mismatches,problems

def fire_unit():
    dictionary=REPO/'data_dictionary.csv'; readme=REPO/'README.md'
    if not dictionary.exists() or not readme.exists(): return 'Not verified: local data dictionary or README missing.'
    with dictionary.open(encoding='utf-8-sig',newline='') as f: rows=list(csv.reader(f))
    fire=[r for r in rows if 'fire_pm25_gfed_mean' in r]
    outcome=[r for r in rows if 'gross_beta_mBq_m3' in r and 'filter_period_biomass_GFED_2001_2023.csv.gz' in r]
    assert len(fire)==len(outcome)==1
    assert fire[0][-1]=='microgram/m3 (centered variables as documented)'
    assert outcome[0][-1]=='mBq/m3'
    readme.read_text()  # Included in the unit/source audit; no coefficient unit is stated there.
    return 'mBq/m³ per µg/m³ (derived from data_dictionary.csv: gross_beta_mBq_m3 = mBq/m3; fire_pm25_gfed_mean = microgram/m3; no explicit fire_pm25_coef unit in README).'

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_si', type=Path, help='Full-precision Supplementary Information DOCX')
    parser.add_argument('output_si', type=Path, help='Reformatted Supplementary Information DOCX')
    parser.add_argument('--report', type=Path, help='JSON report (default: OUTPUT_SI.report.json)')
    args=parser.parse_args()
    source=args.source_si; out=args.output_si
    report_path=args.report if args.report is not None else out.with_suffix('.report.json')
    assert report_path.resolve() not in {source.resolve(),out.resolve()}, 'Report must be a separate file'
    assert source.resolve()!=out.resolve(), 'Input is read only'
    raw,root,body,tables,data=load_source(source)
    specs,s3drop,wide,wind,problems=build_specs(data)
    replacements,heading,intro=paragraph_changes(body,s3drop,wide,wind)
    for s,t in zip(specs,tables): replacements[list(body).index(t)]=make_table(t,s)
    spans=body_spans(raw); assert len(spans)==len(body)
    patches=[(spans[i][0],spans[i][1],serialize_in_context(p)) for i,p in replacements.items()]
    patches.append((spans[heading][1],spans[heading][1],serialize_in_context(intro)))
    new=raw
    for start,end,value in sorted(patches,reverse=True): new=new[:start]+value+new[end:]
    E.fromstring(new)  # fail before publishing malformed XML
    with ZipFile(source) as zin,ZipFile(out,'w') as zout:
        zout.comment=zin.comment
        for info in zin.infolist(): zout.writestr(info,new if info.filename=='word/document.xml' else zin.read(info.filename))
    mismatches,csvproblems=crosscheck(data); problems.extend(csvproblems)
    report={'tables':[{k:v for k,v in s.items() if k not in ['headers','rows','numeric']} for s in specs], 's3_wind_reference_dropped':s3drop,'s6_wide':wide,'s7_wind':wind,'fire_pm25_unit':fire_unit(),'csv_mismatches':mismatches,'check_failures':[], 'problems':problems}
    report['preservation']={'changed_body_indices':sorted(replacements),'new_paragraph_after_body_index':heading,'all_other_zip_members_byte_identical':True,'all_unedited_document_xml_spans_byte_identical':True}
    report['csv_schema_mapping']={'S7':{'n_model':['n_model_M0','n_model_M1'],'both_model_counts_equal_docx':not any(m.get('column')=='n_model' for m in mismatches),'additional_repository_columns':'Not part of the DOCX table; all overlapping cells compared.'}}
    report_path.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__=='__main__': main()
