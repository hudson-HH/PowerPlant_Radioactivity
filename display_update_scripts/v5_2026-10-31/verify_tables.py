"""Compare printed SI numbers with the original full-precision SI DOCX.

Usage: python3 verify_tables.py SOURCE_SI.docx PRINTED_SI.docx
"""
import argparse
import zipfile, re, math
from decimal import Decimal
from lxml import etree
W='{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
def tabs(f):
    d=etree.fromstring(zipfile.ZipFile(f).read('word/document.xml'))
    out=[]
    for t in d.find(W+'body').iter(W+'tbl'):
        out.append([[''.join(x.text or '' for x in c.iter(W+'t')) for c in r.findall(W+'tc')] for r in t.findall(W+'tr')])
    return out
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('source_si', help='Original full-precision Supplementary Information DOCX')
parser.add_argument('printed_si', help='Reformatted Supplementary Information DOCX')
args=parser.parse_args()
A=tabs(args.source_si); B=tabs(args.printed_si)
NUM=re.compile(r'<?−?[\d,]*\.?\d+')
def nums(s): return [float(x.replace('−','-').replace(',','').lstrip('<')) for x in NUM.findall(s)]
def sig3(x):
    if x==0: return 0.0
    return float(f'{x:.3g}')
def ok(src, shown, kind):
    if kind=='int': return int(float(src))==shown
    if kind=='w': return round(float(src),3)==shown
    if kind=='p': return (float(src)<0.001) if shown==0.001 else round(float(src),3)==shown
    return sig3(float(src))==shown
bad=0; n=0
def chk(ti, row_src, row_new, spec):
    """spec: list of (new_col, [src cols], kind)"""
    global bad, n
    for nc, scs, kind in spec:
        vals=nums(row_new[nc])
        if '<' in row_new[nc]: vals=[0.001]
        for sc, v in zip(scs, vals):
            n+=1
            if not ok(row_src[sc], v, kind): bad+=1; print('MISMATCH', ti, row_src[:4], row_new[nc], row_src[sc])
        if len(vals)!=len(scs): bad+=1; print('COUNT', ti, row_new[nc], scs)
    if any('-' in c for c in row_new if re.search(r'\d',c)): print('ASCII hyphen', ti, row_new)
H=lambda t: {h:i for i,h in enumerate(t[0])}
specs={
0:lambda h:[(3,[h['n_rows']],'int'),(4,[h['n_monitors']],'int'),(5,[h['n_exposed_rows']],'int'),(6,[h['n_exposed_monitors']],'int'),(7,[h['n_monitor_years_with_variation']],'int')],
1:lambda h:[(2,[h['n_rows']],'int'),(3,[h['n_monitors']],'int'),(4,[h['n_exposed_rows']],'int'),(5,[h['n_exposed_monitors']],'int'),(6,[h['n_monitor_years_with_variation']],'int')],
2:lambda h:[(3,[h['Wind_speed_m_s']],'w'),(4,[h['Estimate_mBq_m3'],h['CI95_lower'],h['CI95_upper']],'s'),(5,[h['N_pre']],'int'),(6,[h['N_model']],'int')],
3:lambda h:[(3,[h['Wind_speed_m_s']],'w'),(4,[h['N_monitors']],'int'),(5,[h['Estimate_mBq_m3'],h['CI95_lower'],h['CI95_upper']],'s'),(6,[h['N_pre']],'int')],
4:lambda h:[(2,[h['Sector_width_deg']],'int'),(3,[h['Estimate_mBq_m3'],h['CI95_lower'],h['CI95_upper']],'s')],
5:lambda h:[(2,[h['estimate'],h['ci95_low'],h['ci95_high']],'s'),(3,[h['se']],'s'),(4,[h['p_wald']],'p')],
7:lambda h:[(2,[h['n_model']],'int'),(3,[h['original_estimate'],h['original_ci95_low'],h['original_ci95_high']],'s'),(4,[h['fire_adjusted_estimate'],h['fire_adjusted_ci95_low'],h['fire_adjusted_ci95_high']],'s'),(5,[h['relative_change_pct']],'s'),(6,[h['fire_pm25_coef'],h['fire_pm25_ci95_low'],h['fire_pm25_ci95_high']],'s'),(7,[h['fire_pm25_p']],'p')],
}
for ti,f in specs.items():
    a,b=A[ti],B[ti]; assert len(a)==len(b),(ti,len(a),len(b))
    for rs,rn in zip(a[1:],b[1:]): chk(ti,rs,rn,f(H(a)))
# S6: pair current/natural_spline
a,b=A[6],B[6]; h=H(a)
cur=[r for r in a[1:] if r[h['model']]=='current']; spl=[r for r in a[1:] if r[h['model']]=='natural_spline']
assert len(cur)==len(spl)==len(b)-1
for c,s,rn in zip(cur,spl,b[1:]):
    assert [c[h[k]] for k in ('fuel','metric','radius_km')]==[s[h[k]] for k in ('fuel','metric','radius_km')]
    chk(6,c,rn,[(3,[h['N_model']],'int'),(4,[h['wind_ms']],'w'),(5,[h['estimate'],h['ci95_low'],h['ci95_high']],'s')])
    chk(6,s,rn,[(6,[h['estimate'],h['ci95_low'],h['ci95_high']],'s')])
    if rn[0].lower().replace('natural gas','gas')!=c[h['fuel']]: print('label', rn[0], c[h['fuel']])
print('checked',n,'numbers; mismatches',bad)
for i,t in enumerate(B): print(i, t[0]); print('   ', t[1])
