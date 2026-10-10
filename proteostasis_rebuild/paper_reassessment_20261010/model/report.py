"""Regenerate RESULTS from saved rows; no hand-entered numerical conclusions."""
import csv
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent

def read(path):
    with path.open() as f:
        return [{k:(v if k in ('model','scenario','control') else float(v)) for k,v in r.items()}
                for r in csv.DictReader(f)]

def report(out):
    rows=read(out/'trials.csv'); controls=read(out/'controls.csv'); sens=read(out/'sensitivity.csv')
    manifest=json.loads((out/'manifest.json').read_text()); n=manifest['trials']
    lookup={(int(r['trial']),r['model'],r['scenario']):r for r in rows}
    get=lambda m,s,key:np.array([lookup[i,m,s][key] for i in range(n)])
    q=lambda a:'/'.join(f'{(0. if abs(x)<1e-12 else x):.6g}' for x in np.quantile(a,[.05,.5,.95]))
    text=['# Results: finite-capacity client routes', '',
          '**Empirical predictive claim: FAIL / not established.** This is a validated computational discrimination of model assumptions, not a fitted E. coli failure probability. GroEL obligation and TF/DnaK route coverage use the observations being reproduced. Neither is held-out validation.', '',
          f'{n} independent seeded parameter draws (seed {manifest["seed"]}); the same draws, demand and total scalar capacities are used across all three nested models. Primary load is 0.9 times the full architecture’s all-client productive radial capacity. This guarantees a feasible baseline by design, not by observation.', '',
          '## Matched single-machine losses', '',
          'Each target loses δ = 0.8 min(C_DnaK,C_GroEL,C_ClpB). Thus remaining scalar capacity is exactly equal and no capacity is created. Values below are fractions of total client demand; ranges are empirical 5th/50th/95th percentiles of the specified generator, not biological confidence intervals.', '',
          '| Model | Baseline folded fraction, median | GroEL loss: folding drop q05/q50/q95 | DnaK loss: folding drop q05/q50/q95 | ClpB loss: folding drop q05/q50/q95 | Nonzero target contrast |',
          '|---|---:|---|---|---|---:|']
    gaps=[]
    for m in ('scalar','independent','architecture'):
        base=get(m,'baseline','folded_fraction')
        drops=[base-get(m,'matched_'+s,'folded_fraction') for s in ('GroEL','DnaK','ClpB')]
        contrast=np.ptp(np.vstack(drops),axis=0)
        text.append(f'| {m} | {np.median(base):.6g} | '+ ' | '.join(q(a) for a in drops)+f' | {np.sum(contrast>1e-7)}/{n} |')
        gaps.append(contrast)
    rem=np.array([lookup[i,'architecture','baseline']['total_capacity']-lookup[i,'architecture','matched_ClpB']['total_capacity'] for i in range(n)])
    text+=['',f'Matched absolute work loss q05/q50/q95: {q(rem)} (baseline total capacity 1.1). The small ClpB pool often limits common perturbation size; a weak contrast here does not identify interchangeable machines.']
    text+=['','## Pair interactions','','For folded fraction F, excess joint loss is I = F(a)+F(b)−F(ab)−F(0), positive for supra-additive loss. Each endpoint loses 0.4 min(Ca,Cb). The ab, aa and bb arms all lose the same total capacity. aa/bb concentrate twice the endpoint dose on one endpoint. These are finite-dose interactions, not genetic-interaction scores.','','| Pair | Model | I q05/q50/q95 | I > 1e-7 | ab/aa/bb contrast > 1e-7 |','|---|---|---|---:|---:|']
    for pair in ('DnaK_GroEL','DnaK_ClpB','GroEL_ClpB','TF_DnaK'):
        for m in ('scalar','independent','architecture'):
            f=lambda arm:get(m,'pair_'+pair+'_'+arm,'folded_fraction')
            interaction=f('a')+f('b')-f('ab')-get(m,'baseline','folded_fraction')
            contrast=np.ptp(np.vstack([f('ab'),f('aa'),f('bb')]),axis=0)
            text.append(f'| {pair} | {m} | {q(interaction)} | {sum(interaction>1e-7)}/{n} | {sum(contrast>1e-7)}/{n} |')
    text+=['','Positive interactions can arise in a scalar model from an ordinary capacity threshold, even when none occur in this primary grid. They cannot establish a biological serial mechanism. Matched arm contrasts can reject the scalar reduction within a specified network, but do not validate that network.','','## Obligation and upstream-overlap ablations','','Complete knockouts below are **not** matched-total perturbations. Demand remains fixed to the original architecture baseline, including in ablations. All bypasses have finite machine costs.','','| Route assumption | Perturbation | Essential fraction q05/q50/q95 | Folding fraction q05/q50/q95 | Essential quota shortfall |','|---|---|---|---|---:|']
    for s in ('KO_GroEL','KO_TF','KO_DnaK','KO_TF_DnaK'):
        e=get('architecture',s,'essential_fraction'); f=get('architecture',s,'folded_fraction')
        text.append(f'| full | {s} | {q(e)} | {q(f)} | {int(sum(get("architecture",s,"essential_shortfall")))}/{n} |')
    for control,ss in [('remove_obligate',['KO_GroEL']),('remove_upstream',['KO_TF','KO_DnaK','KO_TF_DnaK']),('no_serial',['KO_TF_DnaK'])]:
        for s in ss:
            rr=[r for r in controls if r['control']==control and r['scenario']==s]
            text.append(f'| {control} | {s} | {q([r["essential_fraction"] for r in rr])} | {q([r["folded_fraction"] for r in rr])} | {int(sum(r["essential_shortfall"] for r in rr))}/{n} |')
    text+=['','**Tautology audit.** Zero essential productive flux after GroEL knockout is a theorem of excluding every GroEL-free productive route for the essential-obligate class. Likewise, loss of overlap-class rescue after deleting both TF and DnaK follows from specifying those as its only upstream alternatives. Their reproduction is encoded. Single-to-double magnitudes, competition, and which finite-dose perturbation is worst are conditional computations; these are not yet independent empirical predictions. Calloni also reports impaired GroEL function, ribosome biogenesis and aggregation; these effects are not represented by this static LP. Removing upstream dependence only for overlap clients does not remove the separate upstream requirement imposed on obligate clients; no_serial isolates that latter assumption.','','## Load and catalytic-range sensitivity','','All rows use the same 200 trial identities when run with the default configuration; load is recalibrated to the width-specific unperturbed radial capacity. Width w means independent log10 multipliers Uniform[−w,w], not inferred kinetic distributions.','','| Width | Load | Baseline essential shortfalls | Matched target contrast | GroEL worst (strict) | DnaK worst (strict) | ClpB worst (strict) |','|---|---:|---:|---:|---:|---:|']
    for w in (1.,2.):
        for l in (.6,.9,1.2):
            rr={(int(r['trial']),r['scenario']):r for r in sens if r['width']==w and r['load']==l}
            arr=np.array([[rr[i,'matched_'+s]['folded_fraction'] for i in range(n)] for s in ('GroEL','DnaK','ClpB')])
            wins=[sum(arr[j]+1e-7 < np.min(np.delete(arr,j,axis=0),axis=0)) for j in range(3)]
            text.append(f'| {w:g} | {l:g} | {int(sum(rr[i,"baseline"]["essential_shortfall"] for i in range(n)))}/{n} | {sum(np.ptp(arr,axis=0)>1e-7)}/{n} | '+ ' | '.join(f'{v}/{n}' for v in wins)+' |')
    text+=['','## Distinct outcomes and verification','','Essential shortfall means failure to deliver the specified productive quota for the modeled essential class. It is not cell death. Degradation and finite-horizon sequestration can reduce unresolved influx but cannot satisfy that quota. Unserved flux is residual proteotoxic influx, not a measured proteotoxic-failure threshold; no aggregation dynamics or growth law is fitted. ATP is a synthetic route cost, not a measured cellular energy budget. Its minimization is fourth in the lexicographic objective, so it cannot replace productive rescue with cheaper disposal.','','| Architecture arm | Degradation q05/q50/q95 | Sequestration q05/q50/q95 | Unserved q05/q50/q95 | ATP q05/q50/q95 |','|---|---|---|---|---|']
    for s in ('baseline','matched_GroEL','matched_DnaK','matched_ClpB','KO_TF_DnaK'):
        text.append(f'| {s} | '+' | '.join(q(get('architecture',s,k)) for k in ('degradation','sequestration','unserved','atp'))+' |')
    # Exact invariants across the completed ensemble; fail instead of publishing contradictions.
    assert max(gaps[0])<1e-7, 'scalar target identity leakage'
    for i in range(n):
        for s in {r['scenario'] for r in rows}:
            a,b,c=[lookup[i,m,s] for m in ('scalar','independent','architecture')]
            # Essential primary optimum nests. Secondary optimum nests only when
            # essential optima match; don't assert an invalid lexicographic order.
            assert a['essential_fraction']+1e-7>=b['essential_fraction']>=c['essential_fraction']-1e-7
            if abs(a['essential_fraction']-c['essential_fraction'])<1e-7:
                assert a['folded_fraction']+1e-7>=b['folded_fraction']>=c['folded_fraction']-1e-7
        assert np.ptp([lookup[i,'architecture','matched_'+s]['total_capacity'] for s in ('GroEL','DnaK','ClpB')])<1e-12
    for r in rows+controls+sens:
        assert -1e-8 <= r['folded_fraction'] <= 1.+1e-8
        assert -1e-8 <= r['essential_fraction'] <= 1.+1e-8
    maxerr=max(r['certificate_error'] for r in rows+controls+sens)
    text+=['',f'Every route-allocation optimization checks primal feasibility, dual feasibility, nonnegative prices and strong duality. Maximum absolute certificate residual across all saved evaluations: **{maxerr:.3g}**. Ensemble assertions pass: scalar target invariance, equal remaining capacities and nested essential optima (secondary optima checked when primary quotas coincide).', '',
           'The supported mathematical result is an explicit limitation of scalar reductions; the biological numerical ranking is not established. No pass is awarded for matching the observations used to specify routes. Required next evidence: client-resolved synthesis/turnover demands, productive route qualification, machine occupancy/service times, and held-out perturbation outcomes. Kerner client counts and capacity occupancy, Calloni interactions, and Babu genetic interactions are not full rate constraints.', '',
           'Two direct confrontations with existing experiments — the imposed serial TF/DnaK requirement versus GroEL/ES overproduction rescue of the double knockout, and the effect of all-client radial normalization — are in [CONFRONTATION.md](CONFRONTATION.md), produced by `confront.py`.', '',
           'See [THEORY.md](THEORY.md) for proofs, [README.md](README.md) for exact reproduction and scope, and [SOURCES.md](SOURCES.md) for provenance. Raw outputs and input draws are in `outputs/`.']
    target=HERE/'RESULTS.md' if out.resolve()==(HERE/'outputs').resolve() else out/'RESULTS.md'
    target.write_text('\n'.join(text)+'\n')

if __name__=='__main__': report(HERE/'outputs')
