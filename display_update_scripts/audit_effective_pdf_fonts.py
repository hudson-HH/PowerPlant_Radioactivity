#!/usr/bin/env python3
"""Read-only final-PDF transformed font and editable-text gate (no rendering)."""
import argparse
import collections
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import fitz


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    parser.add_argument("svg", type=Path)
    parser.add_argument("--json-out", required=True, type=Path)
    parser.add_argument("--min-pt", default=5.0, type=float)
    args = parser.parse_args()
    sizes, texts, fonts = [], [], set()
    with fitz.open(args.pdf) as doc:
        pages = len(doc)
        dimensions = [[page.rect.width, page.rect.height] for page in doc]
        for page in doc:
            for block in page.get_text("dict")["blocks"]:
                if block.get("type") != 0:
                    continue
                for line in block["lines"]:
                    for span in line["spans"]:
                        if span["text"].strip():
                            sizes.append(span["size"])
                            fonts.add(span["font"])
                            texts.append(span["text"])
    svg = ET.parse(args.svg).getroot()
    svg_texts = [node for node in svg.iter() if node.tag.rsplit("}", 1)[-1] == "text"]
    below = sum(size < args.min_pt - 1e-6 for size in sizes)
    result = {
        "pdf": str(args.pdf), "svg": str(args.svg), "page_count": pages,
        "page_dimensions_pt": dimensions,
        "font_measurement": "PyMuPDF effective span sizes after PDF text-matrix transformation; raw Tf alone is not the effective point size",
        "effective_minimum_pt": min(sizes) if sizes else None,
        "minimum_required_pt": args.min_pt, "below_minimum_count": below,
        "effective_sizes_pt": dict(sorted(collections.Counter(round(s, 4) for s in sizes).items())),
        "font_names": sorted(fonts), "extractable_pdf_text_spans": len(texts),
        "editable_svg_text_elements": len(svg_texts),
        "all_expected_panel_titles_extractable": all(any(token in t for t in texts) for token in ("Common", "Within", "PM2.5")),
        "verdict": "PASS" if sizes and not below and svg_texts else "FAIL",
    }
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False))
    if result["verdict"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
