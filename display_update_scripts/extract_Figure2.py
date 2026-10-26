"""Extract the user's original third slide, without redrawing its contents.

The input PDF is a headless LibreOffice export of the preserved source PPTX.
PyMuPDF preserves vectors/text; Pillow supplies lossless high-resolution TIFF.
"""
from pathlib import Path
import csv, hashlib, json, zipfile
from lxml import etree
import fitz
from PIL import Image

P = Path(__file__).resolve().parents[1]
V = 'figure2'
W = P/'outputs'/'display_updates'/V
D = P/'outputs'/'display_updates'/V
L = P/'validation'/'display_updates'/V
ppt = P/'data'/'external_inputs'/'study_design.pptx'
source_pdf = P/'data'/'external_inputs'/'study_design.pdf'
stem = D/'figures/Figure_2_study_design_framework'
for directory in (W, L, D/'source_data', stem.parent):
    directory.mkdir(parents=True, exist_ok=True)

ns={'p':'http://schemas.openxmlformats.org/presentationml/2006/main',
    'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
    'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
with zipfile.ZipFile(ppt) as z:
    deck=etree.fromstring(z.read('ppt/presentation.xml'))
    relmap={r.get('Id'):r.get('Target') for r in etree.fromstring(z.read('ppt/_rels/presentation.xml.rels'))}
    ids=deck.xpath('//p:sldId/@r:id',namespaces=ns)
    member='ppt/'+relmap[ids[2]]
    xml=z.read(member)
    root=etree.fromstring(xml)
    strings=root.xpath('//a:t/text()',namespaces=ns)
    shapes=root.findall('.//p:sp',ns)
    pictures=root.findall('.//p:pic',ns)
    (D/'source_data/Figure_2_original_slide3.xml').write_bytes(xml)
    shape_rows=[]
    for s in shapes:
        props=s.find('p:nvSpPr/p:cNvPr',ns)
        shape_rows.append({'shape_id':props.get('id'),'shape_name':props.get('name'),
                           'text':' | '.join(s.xpath('.//a:t/text()',namespaces=ns))})
with (D/'source_data/Figure_2_original_slide3_shape_inventory.csv').open('w',newline='') as f:
    out=csv.DictWriter(f,fieldnames=['shape_id','shape_name','text'])
    out.writeheader();out.writerows(shape_rows)

pdf=fitz.open(source_pdf)
assert len(pdf)==len(ids), 'Unexpected source PDF slide count'
original=fitz.open();original.insert_pdf(pdf,from_page=2,to_page=2)
original.save(W/'Figure_2_original_slide3.pdf')
# Uniform scaling only; no cropping, relabeling or content modifications.
width=183/25.4*72
height=width*pdf[2].rect.height/pdf[2].rect.width
final=fitz.open();page=final.new_page(width=width,height=height)
page.show_pdf_page(page.rect,pdf,2)
final.save(stem.with_suffix('.pdf'))
page.get_pixmap(matrix=fitz.Matrix(600/72,600/72),alpha=False).save(stem.with_suffix('.png'))
with Image.open(stem.with_suffix('.png')) as im:
    im.save(stem.with_suffix('.tiff'),compression='tiff_lzw',dpi=(600,600))
stem.with_suffix('.svg').write_text(page.get_svg_image(text_as_path=False))
normalize=lambda s:''.join(s.split()).replace('–','-').replace('—','-')
all_pdf=normalize(page.get_text())
missing=[s for s in strings if normalize(s) not in all_pdf]
assert not missing, missing
spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for line in b['lines'] for s in line['spans']]
proof={'original_pptx':str(ppt),'pptx_sha256':hashlib.sha256(ppt.read_bytes()).hexdigest(),
       'slide_number':3,'member':member,'slides':len(ids),'native_shapes':len(shapes),'pictures':len(pictures),
       'operation':'Direct vector PDF export; uniform scaling to 183 mm width; no content edits or crop',
       'width_mm':width/72*25.4,'height_mm':height/72*25.4,'dpi':600,
       'text_strings_expected':len(strings),'missing_text':missing,
       'effective_font_pt':sorted(set(round(s['size'],3) for s in spans)),
       'sources':strings,'output_pdf_sha256':hashlib.sha256(stem.with_suffix('.pdf').read_bytes()).hexdigest()}
(L/'Figure_2_extraction_validation.json').write_text(json.dumps(proof,ensure_ascii=False,indent=2))
print(json.dumps({k:v for k,v in proof.items() if k!='sources'},ensure_ascii=False,indent=2))
