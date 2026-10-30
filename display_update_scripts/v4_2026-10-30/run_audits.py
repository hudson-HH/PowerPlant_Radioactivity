"""Inspect final PDF/raster outputs and run the backend-neutral geometry tools."""
from pathlib import Path
import json,subprocess,sys,hashlib
import fitz
from PIL import Image
ROOT=Path(__file__).resolve().parents[2]
summ={}
for pdf in sorted((ROOT/'data/display_assets/v4_2026-10-30/figures').glob('*.pdf')):
    fid=next(k for k,v in json.loads((ROOT/'display_update_scripts/v4_2026-10-30/figure_names.json').read_text()).items() if v==pdf.stem)
    doc=fitz.open(pdf);page=doc[0]
    font_rows=[]
    for font in page.get_fonts(full=True):
        name,ext,typ,data=doc.extract_font(font[0])
        font_rows.append(dict(font=font[3],embedded=bool(data),bytes=len(data)))
    all_embedded=all(f['embedded'] for f in font_rows)
    assert all_embedded,(fid,font_rows)
    png=Image.open(pdf.with_suffix('.png'));tiff=Image.open(pdf.with_suffix('.tiff'))
    assert png.size==tiff.size
    assert all(abs(x-600)<.1 for x in png.info['dpi'])
    assert all(abs(x-600)<.1 for x in tiff.info['dpi'])
    assert tiff.info['compression']=='tiff_lzw'
    # All R scientific plots must remain vector: no raster XObject.
    if fid!='Figure_2':assert len(page.get_images())==0
    row=dict(pdf_size_pt=list(page.rect),png_pixels=png.size,png_dpi=png.info['dpi'],
      tiff_dpi=[float(x) for x in tiff.info['dpi']],tiff_compression=tiff.info['compression'],
      vector_scientific_layers=(fid!='Figure_2'),images=len(page.get_images()),all_fonts_embedded=all_embedded,fonts=font_rows,text=page.get_text())
    for tool,args in [
      ('audit_figure_collisions.py',[str(pdf),'--json-out',str(ROOT/'data/display_assets/v4_2026-10-30/qa/audits'/f'{fid}.collision.json'),'--overlay-pdf',str(ROOT/'data/display_assets/v4_2026-10-30/qa/audits'/f'{fid}.collision.pdf')]),
      ('audit_pdf_text.py',[str(pdf),'--min-pt','5','--json']),
    ]:
        p=subprocess.run([sys.executable,str(ROOT/'display_update_scripts/v4_2026-10-30/qa_helpers'/tool),*args],capture_output=True,text=True)
        tag='collision' if 'collision' in tool else 'fonts'
        (ROOT/'data/display_assets/v4_2026-10-30/qa/logs'/f'{fid}.{tag}.log').write_text(p.stdout+p.stderr)
        row[tag+'_exit']=p.returncode
    layout=ROOT/'data/display_assets/v4_2026-10-30/qa/audits'/f'{fid}.alignment-layout.json'
    if layout.exists():
        p=subprocess.run([sys.executable,str(ROOT/'display_update_scripts/v4_2026-10-30/qa_helpers/audit_panel_alignment.py'),str(layout),
          '--json-out',str(ROOT/'data/display_assets/v4_2026-10-30/qa/audits'/f'{fid}.alignment.json'),'--overlay-svg',str(ROOT/'data/display_assets/v4_2026-10-30/qa/audits'/f'{fid}.alignment.svg'),'--strict'],capture_output=True,text=True)
        (ROOT/'data/display_assets/v4_2026-10-30/qa/logs'/f'{fid}.alignment.log').write_text(p.stdout+p.stderr)
        row['alignment_exit']=p.returncode
    summ[fid]=row
    print(fid,'fonts embedded',all_embedded,'collisions',row['collision_exit'],'alignment',row.get('alignment_exit'))
    (ROOT/'data/display_assets/v4_2026-10-30/qa/audits/export_validation.json').write_text(json.dumps(summ,indent=2))
