"""Run from any directory. All artifacts stay beside this script."""
import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path
import time
import numpy as np
import scipy
from model import architecture, draw, evaluate, radial_capacity, perturbations, pools, MACHINES, ROOT

HERE=Path(__file__).resolve().parent
KINDS=('scalar','independent','architecture')

def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def run(n=200, seed=20261010, out=HERE/'outputs'):
    out.mkdir(parents=True,exist_ok=True)
    rows=[]; controls=[]; inputs=[]; sensitivity=[]
    started=time.monotonic()
    for i in range(n):
        shares,C,costs,multipliers=draw(seed,i)
        net=architecture(costs); lam=radial_capacity(net,shares,C); d=.9*lam*shares
        inputs.append(dict(trial=i, shares=shares.tolist(), capacity=C.tolist(), costs=costs.tolist(),
                           catalytic_multipliers=multipliers.tolist(), radial_capacity=lam))
        for scenario,c in perturbations(C).items():
            for kind in KINDS:
                metrics=evaluate(net,d,c,kind)
                rows.append(dict(trial=i, model=kind, scenario=scenario, total_capacity=float(c.sum()),
                                 demand=float(d.sum()), **metrics))
        for control, scenarios in [('remove_obligate',('KO_GroEL',)),
                                    ('remove_upstream',('KO_TF','KO_DnaK','KO_TF_DnaK')),
                                    ('no_serial',('KO_TF_DnaK',))]:
            alt=architecture(costs,control)
            for scenario in ('baseline',)+scenarios:
                metrics=evaluate(alt,d,perturbations(C)[scenario],'architecture')
                controls.append(dict(trial=i,control=control,scenario=scenario,**metrics))
        # Same random draws across width/load settings. Reference demand is
        # recalibrated per width, NEVER per perturbed machine or model.
        for width in (1.,2.):
            s,c,cc,_=draw(seed,i,width); nn=architecture(cc); ll=radial_capacity(nn,s,c)
            for load in (.6,.9,1.2):
                for scenario,cp in list(perturbations(c).items())[:4]:
                    metrics=evaluate(nn,load*ll*s,cp,'architecture')
                    sensitivity.append(dict(trial=i,width=width,load=load,scenario=scenario,**metrics))
        if (i+1)%25==0:
            print(f'{i+1}/{n} trials, {time.monotonic()-started:.1f}s', flush=True)
    write_csv(out/'trials.csv', rows); write_csv(out/'controls.csv',controls)
    write_csv(out/'sensitivity.csv',sensitivity)
    (out/'inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')
    manifest=dict(seed=seed,trials=n,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                  pool_units='uM functional units; relative scale only',pools=dict(zip(MACHINES[:7],pools())),
                  hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in [ROOT/'WO-07/bundles.tsv',HERE/'model.py',HERE/'run.py',HERE/'report.py']})
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    from report import report
    report(out)

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--trials',type=int,default=200)
    p.add_argument('--seed',type=int,default=20261010); p.add_argument('--out',type=Path,default=HERE/'outputs')
    a=p.parse_args(); run(a.trials,a.seed,a.out)
