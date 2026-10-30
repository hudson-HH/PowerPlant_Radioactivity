"""Read-only font audit using transformed PDF glyph sizes, not raw Tf operands."""
from pathlib import Path
import json,collections
import fitz
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[3]
rows={}
for pdf in sorted((ROOT/'data/display_assets/v4_2026-10-30/figures').glob('*.pdf')):
    sizes=[];names=set();text=[]
    doc=fitz.open(pdf)
    for page in doc:
        for block in page.get_text('dict')['blocks']:
            for line in block.get('lines',[]):
                for span in line['spans']:
                    if span['text'].strip():
                        sizes.append(span['size']);names.add(span['font']);text.append(span['text'])
    svg=ET.parse(pdf.with_suffix('.svg')).getroot()
    editable=sum(n.tag.rsplit('}',1)[-1]=='text' for n in svg.iter())
    rows[pdf.stem]=dict(measurement='Effective PDF span sizes after text-matrix transforms',
      effective_minimum_pt=min(sizes),sizes_pt=dict(collections.Counter(round(x,4) for x in sizes)),
      font_names=sorted(names),editable_svg_text_elements=editable,
      verdict='PASS' if sizes and min(sizes)>=5-1e-6 and editable else 'FAIL')
    print(pdf.stem,rows[pdf.stem]['verdict'],round(min(sizes),3),'pt')
(ROOT/'data/display_assets/v4_2026-10-30/qa/audits/effective_fonts.json').write_text(json.dumps(rows,indent=2))
assert all(r['verdict']=='PASS' for r in rows.values())
