from pathlib import Path
import csv, subprocess, sys, tempfile

ROOT=Path(__file__).resolve().parents[2]
CMP=ROOT/"scripts/compare_validation_results.py"
AGG=ROOT/"scripts/aggregate_validation_runs.py"

def write(path, rows):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=["case","Cd","Cl"]); w.writeheader(); w.writerows(rows)

with tempfile.TemporaryDirectory() as d:
    p=Path(d); c=p/"c.csv"; r=p/"r.csv"; o=p/"o.json"
    write(c,[{"case":"naca","Cd":"0.12","Cl":"1e-6"}])
    write(r,[{"case":"naca","Cd":"0.12","Cl":"0"}])
    x=subprocess.run([sys.executable,str(CMP),"--computed",str(c),"--reference",str(r),"--qoi","Cd","Cl","--output",str(o)],capture_output=True,text=True)
    assert x.returncode==0, x.stderr
    assert '"absolute_error": 0.0' in o.read_text()
