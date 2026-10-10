"""AUDIT 2026-10-09 (new work, not part of the original run).

counterfactual architectures for the network-capacity model in ../model.py.
the original solvers are imported unchanged; only the route generator, the
random streams and a few knobs are new.

knobs
  routes     'orig'      the original findRoutesBio rules (reproduced exactly
                         from the same rng draws, see test_audit)
             'density'   same module families, but EVERY optional edge has
                         probability q (q=0: only forced edges; q=1: all edges)
             'nocut'     original rules, then a backup route is added to every
                         class whose routes share a machine, built only from
                         modules disjoint from that shared set. removes every
                         structural single cut by construction
             'unstruct'  original adversarial null (1-4 random routes of 1-3
                         machines)
  agg_free   True: the aggregate class always gets one ClpB-free route
             (SEQ or LON, chosen at random). removes ClpB/KJE exclusivity only
  structure and capacities use SEPARATE rng streams, so capacities can be
  resampled with the route set held fixed (the original code draws both from
  one stream, so every W setting re-uses the same 1000 route sets)
  abundance  'measured' | 'blind' (a=1) | 'shuffled' (measured values permuted
             across machines with the kinetic stream)
  agg_share  None: Dirichlet demand as original. a float f fixes d_AGG = f and
             rescales the other five classes to 1 - f
  cost       'orig' (primary cheapest) | 'flat' (all routes equal, tiny jitter)
             | 'reverse' (primary most expensive). affects only the order in
             which machines saturate below lambda*, never lambda* itself
"""
import itertools
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import model as M  # noqa: E402

FOLD_MODS = ["TF", "KJE", "GRO", "HTPG"]


def machinesOf(mods):
    return set(m for mod in mods for m in M.MODULES[mod])


def findRoutesDensity(rng, q):
    """original module families, every optional edge at probability q.
    forced edges (kept at q=0): NASC<-TF, DK<-KJE, GRO<-GRO, HTP<-HTPG,
    DMG<-one protease, AGG<-CLPB+KJE. these are the minimum the original
    rules always contain."""
    p = lambda: rng.random() < q
    R = [("NASC", ("TF",))]
    R += [("NASC", (mod,)) for mod in ("KJE", "GRO") if p()]
    R.append(("DK", ("KJE",)))
    R += [("DK", (mod,)) for mod in ("GRO", "TF") if p()]
    R.append(("GRO", ("GRO",)))
    R += [("GRO", (mod,)) for mod in ("KJE", "TF") if p()]
    R.append(("HTP", ("HTPG", "KJE") if p() else ("HTPG",)))
    R += [("HTP", (mod,)) for mod in ("KJE", "GRO") if p()]
    first = M.PROTEASES[rng.integers(3)]
    R.append(("DMG", (first,)))
    R += [("DMG", (q2,)) for q2 in M.PROTEASES if q2 != first and p()]
    R.append(("AGG", ("CLPB", "KJE")))
    R += [("AGG", (mod,)) for mod in ("SEQ", "LON") if p()]
    R += [("AGG", ("CLPB", "KJE", q2)) for q2 in M.PROTEASES if p()]
    return R


def addBackups(rng, R):
    """for each class whose routes all share >= 1 machine, add one route made
    of a class-appropriate module that avoids every shared machine."""
    fam = {k: FOLD_MODS for k in M.FOLDING}
    fam["DMG"] = list(M.PROTEASES)
    fam["AGG"] = ["SEQ", "LON"] + list(M.PROTEASES)
    out = list(R)
    for k in M.CLASSES:
        sets = [machinesOf(mods) for kk, mods in out if kk == k]
        common = set.intersection(*sets)
        if common:
            ok = [mod for mod in fam[k] if not (machinesOf((mod,)) & common)]
            out.append((k, (ok[rng.integers(len(ok))],)))
    return out


def findRoutes(rng, routes="orig", q=None, agg_free=False, rx=None):
    """returns (base routes, extra routes). base routes come from the structure
    stream exactly as in the original; extras come from a separate stream rx so
    that U and d of the base routes stay paired with the original instance."""
    if routes in ("orig", "nocut"):
        R = M.findRoutesBio(rng)
    elif routes == "density":
        R = findRoutesDensity(rng, q)
    elif routes == "unstruct":
        R = M.findRoutesUnstructured(rng)
    else:
        raise ValueError(routes)
    X = []
    if agg_free and all("ClpB" in machinesOf(mods) for k, mods in R if k == "AGG"):
        X.append(("AGG", (("SEQ", "LON")[rx.integers(2)],)))
    if routes == "nocut":
        X += addBackups(rx, R + X)[len(R + X):]
    return R, X


def makeInstanceCF(seed, routes="orig", q=None, agg_free=False, W=2.0, wU=1.0,
                   alpha=1.0, abundance="measured", agg_share=None, cost="orig",
                   pools=None, kin_seed=None):
    """structure stream: seed (routes, U, d: identical to model.makeInstance
    for routes='orig'). extra-route stream: seed + 5_000_000. kinetic stream
    (s, abundance shuffle): (kin_seed or seed) + 7_919_000."""
    rs = np.random.default_rng(seed)
    rx = np.random.default_rng(seed + 5_000_000)
    rk = np.random.default_rng((seed if kin_seed is None else kin_seed) + 7_919_000)
    base, extra = findRoutes(rs, routes, q, agg_free, rx)
    R = base + extra
    nR, nM = len(R), len(M.MACHINES)
    U = np.zeros((nR, nM))
    cls = np.zeros(nR, dtype=int)
    for i, (k, mods) in enumerate(R):
        cls[i] = M.CLASSES.index(k)
        g = rs if i < len(base) else rx
        for mod in mods:
            for m in (M.MODULES[mod] if not mod.startswith("@") else (mod[1:],)):
                U[i, M.MIDX[m]] = 10 ** g.uniform(-wU / 2, wU / 2)
        if i == len(base) - 1:
            d = rs.dirichlet(alpha * np.ones(len(M.CLASSES)))
    if agg_share is not None:
        a_ = M.CLASSES.index("AGG")
        rest = np.delete(d, a_)
        d = np.insert(rest / rest.sum() * (1 - agg_share), a_, agg_share)
    a = np.array([pools[m] for m in M.MACHINES])
    if abundance == "blind":
        a = np.ones(nM)
    elif abundance == "shuffled":
        a = rk.permutation(a)
    s = 10 ** rk.uniform(-W / 2, W / 2, size=nM)
    primary = np.array([int(np.flatnonzero(cls == k)[0]) for k in range(len(M.CLASSES))])
    rank = np.zeros(nR)
    for k in range(len(M.CLASSES)):
        idx = np.flatnonzero(cls == k)
        r_ = np.arange(len(idx), dtype=float)
        rank[idx] = r_ if cost == "orig" else (0 * r_ if cost == "flat" else r_[::-1])
    jit = 1e-3 * rx.random(nR)
    return {"seed": seed, "routes": R, "U": U, "cls": cls, "d": d, "c": a * s, "a": a,
            "s": s, "primary": primary, "cost": 1.0 + rank + jit}


def findPairSplit(inst, rho_list=(0.25, 0.5, 0.9)):
    """architecture C, full knockout. splits synthetic-lethal pairs into
    STRUCTURAL (pair knockout leaves some class with no usable route: W2 = 0,
    true at any load) and LOAD-DEPENDENT (0 < W2 < rho)."""
    lam0, W1, W2 = M.findKnockdowns(inst, "C", 0.0)
    out = {"lamC": lam0, "single_cuts": int(np.sum(W1 == 0))}
    for rho in rho_list:
        st = ld = 0
        for (i, j), w in W2.items():
            if W1[i] >= rho and W1[j] >= rho and w < rho - 1e-9:
                if w == 0:
                    st += 1
                else:
                    ld += 1
        out[f"SL_struct_rho{rho}"] = st
        out[f"SL_load_rho{rho}"] = ld
    return out


def findBottleneck(inst):
    """machines with a positive dual price at lambda*_C."""
    _, _, y = M.solveLamC(inst, duals=True)
    yn = y / y.max() if y.max() > 0 else y
    return ";".join(M.MACHINES[i] for i in np.flatnonzero(yn > 1e-6))


def evalInstance(inst, saturation=True):
    row = findPairSplit(inst)
    row["bottlenecks"] = findBottleneck(inst)
    if saturation:
        ev = M.findSaturationEvents(inst, n=100)
        row["sat_events"] = ev["n_events_before_failure"]
        row["sat_order"] = ";".join(f"{m}@{f:.2f}" for m, f in ev["order"])
    return row


def findCutMachines(R, d=None):
    """structural single cuts from the route list alone (no capacities):
    machines present in every route of some class."""
    out = set()
    for k in M.CLASSES:
        sets = [set(m[1:] for m in mods) if mods[0].startswith("@") else machinesOf(mods)
                for kk, mods in R if kk == k]
        if sets:                      # a class absent from R contributes nothing
            out |= set.intersection(*sets)
    return out
