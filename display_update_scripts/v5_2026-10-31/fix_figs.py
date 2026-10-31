"""Tighten figure-to-legend spacing and keep each figure with its caption.

usage: fix_figs.py in.docx out.docx

- Figure paragraphs get single line spacing with no space before/after. Word applies a
  multiple line spacing to the whole image height, so "double" left a gap as tall as the figure.
- A figure paragraph that also holds caption text (main-text Figure 3) is split in two.
- Figure paragraphs styled as Heading 2 lose that style (they showed up in the navigation pane).
- Adds the missing Figure_S1 bookmark used by the SI contents page.
"""
import sys, zipfile, copy, re
from lxml import etree

import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('source_si', help='Input Supplementary Information DOCX (also accepts manuscript DOCX)')
parser.add_argument('output_si', help='Output DOCX')
args=parser.parse_args()
from pathlib import Path
assert Path(args.source_si).resolve()!=Path(args.output_si).resolve(), 'Input is read only'

W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
src, dst = args.source_si, args.output_si
z = zipfile.ZipFile(src)
doc = etree.fromstring(z.read('word/document.xml'))
body = doc.find(W + 'body')
txt = lambda e: ''.join(t.text or '' for t in e.iter(W + 't'))


def ppr(p):
    pr = p.find(W + 'pPr')
    if pr is None:
        pr = etree.SubElement(p, W + 'pPr'); p.remove(pr); p.insert(0, pr)
    return pr


def child(pr, tag, before=()):
    e = pr.find(W + tag)
    if e is None:
        e = etree.Element(W + tag)
        # schema order: pStyle, keepNext, keepLines, ... spacing ... jc ... rPr
        ref = next((pr.find(W + b) for b in before if pr.find(W + b) is not None), None)
        if ref is not None: ref.addprevious(e)
        else:
            rpr = pr.find(W + 'rPr')
            rpr.addprevious(e) if rpr is not None else pr.append(e)
    return e


AFTER_SPACING = ('ind', 'contextualSpacing', 'jc', 'textAlignment', 'outlineLvl', 'rPr')
log = []
paras = [p for p in body if p.tag == W + 'p']
for p in paras:
    if p.find('.//' + W + 'drawing') is None:
        continue
    # split off caption text that shares the figure paragraph
    runs = [r for r in p if r.tag == W + 'r']
    di = next(i for i, r in enumerate(runs) if r.find('.//' + W + 'drawing') is not None)
    tail = [e for e in list(p) if e.tag != W + 'pPr' and list(p).index(e) > list(p).index(runs[di])]
    if txt(p).strip():
        cap_src = next((q for q in paras if re.match(r'Figure \d+\.', txt(q)) and q.find('.//' + W + 'drawing') is None), None)
        newp = etree.Element(W + 'p')
        newp.append(copy.deepcopy(cap_src.find(W + 'pPr')))
        for e in tail: newp.append(e)
        # the text used to inherit bold and 12 pt from Heading 2; write them on the runs
        for r in newp.iter(W + 'r'):
            rp = r.find(W + 'rPr')
            if rp is None: rp = etree.Element(W + 'rPr'); r.insert(0, rp)
            if rp.find(W + 'b') is None:
                b = etree.Element(W + 'b'); bc = etree.Element(W + 'bCs')
                anchor = rp.find(W + 'rFonts')
                (anchor.addnext(b) if anchor is not None else rp.insert(0, b)); b.addnext(bc)
            if rp.find(W + 'sz') is None:
                sz = etree.Element(W + 'sz'); sz.set(W + 'val', '24')
                rp.find(W + 'szCs').addprevious(sz) if rp.find(W + 'szCs') is not None else rp.append(sz)
        p.addnext(newp)
        log.append(f'split caption off figure paragraph: {txt(newp)[:50]!r}')
    pr = ppr(p)
    st = pr.find(W + 'pStyle')
    if st is not None and st.get(W + 'val') == '21':
        pr.remove(st); log.append('removed Heading 2 style from a figure paragraph')
    ol = pr.find(W + 'outlineLvl')
    if ol is not None: pr.remove(ol)
    child(pr, 'keepNext', ('keepLines', 'spacing') + AFTER_SPACING)
    sp = child(pr, 'spacing', AFTER_SPACING)
    for k in list(sp.attrib): del sp.attrib[k]
    sp.set(W + 'before', '0'); sp.set(W + 'after', '0'); sp.set(W + 'line', '240'); sp.set(W + 'lineRule', 'auto')
    jc = child(pr, 'jc', ('textAlignment', 'outlineLvl', 'rPr')); jc.set(W + 'val', 'center')
    log.append(f'figure paragraph set to single spacing')

# empty Heading 2 paragraphs left next to figures
for p in [p for p in body if p.tag == W + 'p']:
    st = p.find(W + 'pPr/' + W + 'pStyle')
    if st is not None and st.get(W + 'val') == '21' and not txt(p).strip() and p.find('.//' + W + 'drawing') is None:
        nb = [q for q in (p.getprevious(), p.getnext()) if q is not None]
        if any(q.find('.//' + W + 'drawing') is not None or re.match(r'Figure \d+\.', txt(q)) for q in nb):
            st.getparent().remove(st); log.append('removed Heading 2 style from an empty paragraph beside a figure')

# SI: figure title above the image should stay on the same page as the image
for p in [p for p in body if p.tag == W + 'p']:
    if re.match(r'Supplementary Fig\. S\d+\.', txt(p)) and p.getnext() is not None and p.getnext().find('.//' + W + 'drawing') is not None:
        child(ppr(p), 'keepNext', ('keepLines', 'spacing') + AFTER_SPACING)

# missing contents-page bookmark for Supplementary Fig. S1
names = {b.get(W + 'name') for b in doc.iter(W + 'bookmarkStart')}
if 'Figure_S1' not in names:
    t = next((p for p in body if p.tag == W + 'p' and txt(p).startswith('Supplementary Fig. S1.')), None)
    if t is not None:
        nid = str(max(int(b.get(W + 'id')) for b in doc.iter(W + 'bookmarkStart')) + 1)
        bs = etree.Element(W + 'bookmarkStart'); bs.set(W + 'id', nid); bs.set(W + 'name', 'Figure_S1')
        be = etree.Element(W + 'bookmarkEnd'); be.set(W + 'id', nid)
        ppr(t).addnext(bs); bs.addnext(be)
        log.append('added Figure_S1 bookmark')

out = zipfile.ZipFile(dst, 'w', zipfile.ZIP_DEFLATED)
for it in z.infolist():
    data = z.read(it.filename)
    if it.filename == 'word/document.xml':
        data = etree.tostring(doc, xml_declaration=True, encoding='UTF-8', standalone=True)
    out.writestr(it, data)
out.close()
from collections import Counter
for k, v in Counter(log).items(): print(v, 'x', k)
