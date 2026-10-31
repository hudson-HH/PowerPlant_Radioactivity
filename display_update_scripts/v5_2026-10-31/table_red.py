"""In tables: header row red, body rows black. Outside tables untouched. usage: in out"""
import sys, zipfile
from lxml import etree
import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('source_si', help='Input Supplementary Information DOCX (also accepts manuscript DOCX)')
parser.add_argument('output_si', help='Output DOCX')
args=parser.parse_args()
from pathlib import Path
assert Path(args.source_si).resolve()!=Path(args.output_si).resolve(), 'Input is read only'

W='{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
z=zipfile.ZipFile(args.source_si); d=etree.fromstring(z.read('word/document.xml'))
def col(r,v):
    rp=r.find(W+'rPr')
    if rp is None: rp=etree.Element(W+'rPr'); r.insert(0,rp)
    c=rp.find(W+'color')
    if c is None:
        c=etree.Element(W+'color'); sz=rp.find(W+'sz')
        (sz.addprevious(c) if sz is not None else rp.append(c))
    c.attrib.clear(); c.set(W+'val',v)
for t in d.iter(W+'tbl'):
    for i,tr in enumerate(t.findall(W+'tr')):
        for r in tr.iter(W+'r'): col(r,'FF0000' if i==0 else '000000')
o=zipfile.ZipFile(args.output_si,'w',zipfile.ZIP_DEFLATED)
for it in z.infolist():
    b=z.read(it.filename)
    if it.filename=='word/document.xml': b=etree.tostring(d,xml_declaration=True,encoding='UTF-8',standalone=True)
    o.writestr(it,b)
o.close()
