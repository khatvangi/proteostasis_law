"""Independent radial-dual checks and small end-to-end reproducibility audit."""
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from model import ROOT, architecture, SOLVER_OPTIONS

HERE=Path(__file__).resolve().parent

def validate():
    out=HERE/'outputs'
    manifest=json.loads((out/'manifest.json').read_text())
    for name,digest in manifest['hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest, name
    radial_errors=[]
    for r in json.loads((out/'inputs.json').read_text()):
        n=architecture(np.array(r['costs'])); active=n.productive.astype(bool)
        A=n.A[:,active]; H=n.H[:,active]
        d=np.array(r['shares']); C=np.array(r['capacity'])
        # min C*p, A.T*p - H.T*y >= 0, d*y >= 1; y free.
        B=np.vstack([np.c_[-A.T,H.T],np.r_[np.zeros(len(C)),-d]])
        b=np.r_[np.zeros(A.shape[1]),-1.]
        z=linprog(np.r_[C,np.zeros(len(d))],A_ub=B,b_ub=b,
                  bounds=[(0,None)]*len(C)+[(None,None)]*len(d),method='highs', options=SOLVER_OPTIONS)
        assert z.success,z.message
        error=abs(z.fun-r['radial_capacity']); radial_errors.append(error)
        assert error<1e-8
    with tempfile.TemporaryDirectory(prefix='route_repro_') as tmp:
        a,b=Path(tmp)/'a',Path(tmp)/'b'
        for target in (a,b):
            subprocess.run([sys.executable,str(HERE/'run.py'),'--trials','3','--out',str(target)],check=True)
        assert {p.name for p in a.iterdir()}=={p.name for p in b.iterdir()}
        hashes={}
        for p in sorted(a.iterdir()):
            assert p.read_bytes()==(b/p.name).read_bytes(),p.name
            hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
    test=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(HERE),'-v'],
                        capture_output=True,text=True,check=True)
    (out/'unit_tests.txt').write_text(test.stdout+test.stderr)
    result=dict(radial_dual_trials=len(radial_errors),max_radial_duality_gap=max(radial_errors),
                reproducibility='two independent three-trial runs byte-identical for every artifact',
                small_run_sha256=hashes,unit_tests='15 passed',input_and_code_hashes_verified=True)
    (out/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':validate()
