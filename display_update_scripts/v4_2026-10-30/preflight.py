"""Static check of each plotting module with its actual shared export context."""
from pathlib import Path
import json,subprocess,sys
ROOT=Path(__file__).resolve().parents[2]
modules=['Figure1_display_v3.R','Figure6_display_v3.R','02_reproduce_figures_tables.R',
 '02a_redraw_overall_inventory.R','02b_redraw_supplementary_diagnostics.R',
 'Supplementary_Figure_S1.R','Supplementary_Figure_S2.R','Supplementary_Figure_S6.R','Supplementary_Figure_S7.R']
context=(ROOT/'display_update_scripts/v4_2026-10-30/export_helpers.R').read_text()
(ROOT/'data/display_assets/v4_2026-10-30/qa/audits').mkdir(parents=True,exist_ok=True)
rows={}
for module in modules:
    # Inspection copy only: functions remain in display_update_scripts/v4_2026-10-30/ and are run by run_all.R.
    p=ROOT/'data/display_assets/v4_2026-10-30/qa/audits'/f'{module}.with_export_context.R'
    p.write_text(context+'\n'+(ROOT/'display_update_scripts/v4_2026-10-30'/module).read_text())
    r=subprocess.run([sys.executable,str(ROOT/'display_update_scripts/v4_2026-10-30/qa_helpers/validate_figure.py'),str(p),'--backend','r','--json'],capture_output=True,text=True)
    rows[module]=json.loads(r.stdout)
    print(module,rows[module]['summary']['counts'])
(ROOT/'data/display_assets/v4_2026-10-30/qa/audits/source_preflight_with_context.json').write_text(json.dumps(rows,indent=2))
assert all(r['summary']['counts'].get('FAIL',0)==0 for r in rows.values())
