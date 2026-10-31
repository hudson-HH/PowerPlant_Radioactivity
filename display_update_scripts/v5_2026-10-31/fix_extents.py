"""Make each inline picture's frame (wp:extent) match the picture's own aspect ratio.
usage: fix_extents.py in.docx out.docx"""
import sys, zipfile, re, struct
from lxml import etree
import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('source_si', help='Input Supplementary Information DOCX (also accepts manuscript DOCX)')
parser.add_argument('output_si', help='Output DOCX')
args=parser.parse_args()
from pathlib import Path
assert Path(args.source_si).resolve()!=Path(args.output_si).resolve(), 'Input is read only'

WP='{http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing}'
A='{http://schemas.openxmlformats.org/drawingml/2006/main}'
R='{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'
z=zipfile.ZipFile(args.source_si); d=etree.fromstring(z.read('word/document.xml'))
rels=z.read('word/_rels/document.xml.rels').decode()
rid={m.group(1):m.group(2) for m in re.finditer(r'Id="([^"]+)"[^>]*Target="([^"]+)"',rels)}
rid.update({m.group(2):m.group(1) for m in re.finditer(r'Target="([^"]+)"[^>]*Id="([^"]+)"',rels)})
for k in d.iter(WP+'inline'):
    blip=k.find('.//'+A+'blip')
    if blip is None: continue
    b=z.read('word/'+rid[blip.get(R+'embed')])
    if b[:8]!=b'\x89PNG\r\n\x1a\n': continue
    w,h=struct.unpack('>II',b[16:24])
    e=k.find(WP+'extent'); cx=int(e.get('cx')); cy=round(cx*h/w)
    old=e.get('cy'); e.set('cy',str(cy))
    for x in k.iter(A+'ext'):
        if x.get('cx'): x.set('cx',str(cx)); x.set('cy',str(cy))
    ee=k.find(WP+'effectExtent')
    if ee is not None:
        for a in 'ltrb': ee.set(a,'0')
    print(f'{rid[blip.get(R+"embed")]}: frame cy {old} -> {cy}')
o=zipfile.ZipFile(args.output_si,'w',zipfile.ZIP_DEFLATED)
for it in z.infolist():
    data=z.read(it.filename)
    if it.filename=='word/document.xml': data=etree.tostring(d,xml_declaration=True,encoding='UTF-8',standalone=True)
    o.writestr(it,data)
o.close()
