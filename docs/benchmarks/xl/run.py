"""run.py TAG ANDES [extra andes args...]: search every dataset of datasets.json found under
$XLBENCH_DATA, score it against its truth, and write $XLBENCH_DATA/results/TAG/summary.json.
ONLY=name1,name2 restricts the datasets. Thermo .raw input needs an andes built with
--features thermo and the .NET 8 runtime (DOTNET_ROOT)."""
import json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import eval as ev

B = os.environ.get('XLBENCH_DATA', '.')
tag, binary, extra = sys.argv[1], sys.argv[2], sys.argv[3:]
only = os.environ['ONLY'].split(',') if os.environ.get('ONLY') else None
os.makedirs(f'{B}/results/{tag}', exist_ok=True)
summary_path = f'{B}/results/{tag}/summary.json'
summary = json.load(open(summary_path)) if os.path.exists(summary_path) else {}
for d in json.load(open(f'{HERE}/datasets.json')):
    raw = f'{B}/raw/{d["raw"]}'
    if not os.path.exists(raw) or (only and d['name'] not in only):
        continue
    out = f'{B}/results/{tag}/{d["name"]}'
    os.makedirs(out, exist_ok=True)
    cmd = [binary, '--spectrum', raw, '--database', f'{B}/{d["fasta"]}', '--xl', d['xl'],
           '--output-pin', f'{out}/{d["name"]}.pin'] + d.get('args', []) + extra
    t = time.time()
    with open(f'{out}/log.txt', 'w') as log:
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        summary[d['name']] = {'error': rc}
        print(d['name'], 'FAILED', rc, flush=True)
        json.dump(summary, open(summary_path, 'w'), indent=1)
        continue
    r = ev.evaluate(f'{out}/{d["name"]}.xl.tsv', d['truth'])
    r['wall_s'] = round(time.time() - t, 1)
    summary[d['name']] = r
    print(d['name'], json.dumps(r), flush=True)
    json.dump(summary, open(summary_path, 'w'), indent=1)
