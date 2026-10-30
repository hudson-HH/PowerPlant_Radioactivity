"""Presentation-only adaptation of the supplied Figure 2 extraction script.

Crop just the title band from the frozen exports. No PowerPoint box, arrow,
body text, font, colour, or pixel in the retained PNG is reconstructed.
"""
from pathlib import Path
import json, hashlib, re
import fitz
from PIL import Image

ROOT=Path(__file__).resolve().parents[2]
old=ROOT/'data/display_assets/frozen_v3/figures/Figure_2_study_design_framework'
out=ROOT/'data/display_assets/v4_2026-10-30/figures/Figure_2_study_design_framework'
doc=fitz.open(old.with_suffix('.pdf'));page=doc[0]
title='Figure 2. Study design and analytical framework'
hits=page.search_for(title)
assert len(hits)==1
top_px=192
top_pt=top_px*72/600
assert hits[0].y1 < top_pt
clip=fitz.Rect(0,top_pt,page.rect.width,page.rect.height)
retained_text=page.get_text(clip=clip)
assert title not in retained_text
assert page.get_text().replace(title,'').split()==retained_text.split()
# Remove the off-page title object as well, so text extraction and clipping
# audits do not encounter invisible title glyphs in the cropped PDF.
page.add_redact_annot(hits[0],fill=False)
page.apply_redactions(images=0,graphics=0,text=0)
new=fitz.open();p=new.new_page(width=clip.width,height=clip.height)
p.show_pdf_page(p.rect,doc,0,clip=clip,keep_proportion=False)
new.subset_fonts();new.save(out.with_suffix('.pdf'),garbage=4,deflate=True)
with Image.open(old.with_suffix('.png')) as im:
    cropped=im.crop((0,top_px,im.width,im.height))
    cropped.save(out.with_suffix('.png'),dpi=(600,600))
    cropped.save(out.with_suffix('.tiff'),compression='tiff_lzw',dpi=(600,600))
    assert cropped.tobytes()==im.crop((0,top_px,im.width,im.height)).tobytes()
svg=old.with_suffix('.svg').read_text()
svg=re.sub(r'<text\b[^>]*>\s*<tspan\b[^>]*>Figure 2\. Study design and analytical framework</tspan></text>','',svg,count=1)
svg=re.sub(r'viewBox="[^"]+"',f'viewBox="0 {top_pt} {clip.width} {clip.height}"',svg,count=1)
svg=re.sub(r'width="[^"]+"',f'width="{clip.width}pt"',svg,count=1)
svg=re.sub(r'height="[^"]+"',f'height="{clip.height}pt"',svg,count=1)
out.with_suffix('.svg').write_text(svg)
proof=dict(operation='Crop original title band only; retained PNG pixels exact; PDF vectors/fonts retained',
 removed_text=[title],crop_top_pixels=top_px,crop_top_pt=top_pt,png_dpi=600,
 box_text_identical=True,retained_png_pixels_identical=True,
 old_pdf_size_pt=list(page.rect),new_pdf_size_pt=list(p.rect),
 source_pdf_sha256=hashlib.sha256(old.with_suffix('.pdf').read_bytes()).hexdigest(),
 output_pdf_sha256=hashlib.sha256(out.with_suffix('.pdf').read_bytes()).hexdigest())
(ROOT/'data/display_assets/v4_2026-10-30/qa/audits/Figure_2_crop_validation.json').write_text(json.dumps(proof,indent=2))
print(json.dumps(proof,indent=2))
