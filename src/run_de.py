"""Resumable multi-seed Differential Evolution scan for the 2HDM+S paper.

Examples:
  python src/run_de.py --preset pilot --targets 30
  python src/run_de.py --preset production --targets 1000
  python src/run_de.py --preset production --targets 1000 --resume
No legacy output is overwritten; each configuration has its own data directory.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import platform
import sys
from pathlib import Path
from datetime import datetime, timezone
import inspect
import numpy as np
import scipy
from scipy.optimize import NonlinearConstraint, differential_evolution
from model import BOUNDS, MASS_CENTERS, MASS_WIDTHS, V_EW, MC_TARGET, KAPPA_TARGET, spectral_loss, constraints, diagnostics, PARAMETERS
try:
    from fast_model import objective_fast, constraint_fast
except ImportError:
    objective_fast, constraint_fast = spectral_loss, constraints

COLUMNS=['target_index','restart','seed','seed_targets','loss','optimizer_success','status','nfev','nit','mass_pass','aux_pass','all_pass','optimizer_constraints_ok',
         *PARAMETERS, 'target_H1','target_H2','target_H3', 'mass_H1','mass_H2','mass_H3','mc','kappa2V','z_H1','z_H2','z_H3','z_charged','z_kappa','minimum_eigenvalue','message']

def atomic_write_rows(path,rows):
    temp=path.with_suffix('.tmp')
    with temp.open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=COLUMNS);w.writeheader();w.writerows(rows)
    temp.replace(path)

def de_random_kwargs(seed):
    """Pass the same numpy Generator to SciPy DE across supported versions.

    SciPy < 1.15 uses ``seed``; newer SciPy uses ``rng``. Older SciPy
    versions also accept numpy.random.Generator as ``seed``, so using a
    Generator in both cases avoids switching to legacy RandomState.
    """
    generator = np.random.default_rng(seed)
    parameter = "rng" if "rng" in inspect.signature(differential_evolution).parameters else "seed"
    return {parameter: generator}


def run_one(target, target_idx, repeat, seed, options):
    # Both constraints are optimized simultaneously: original charged mass^2
    # residual +/- 260 GeV^2; canonical sign-independent kappa in [0.98,1.08].
    nonlinear=NonlinearConstraint(constraint_fast,[-260.0,0.98],[260.0,1.08])
    result=differential_evolution(
        objective_fast, BOUNDS, args=(target,), constraints=(nonlinear,),
        maxiter=options['maxiter'],popsize=options['popsize'],
        tol=options['tol'], mutation=(0.5,1.0),recombination=0.7,
        strategy='best1bin', polish=options['polish'],
        workers=1, updating='immediate', disp=False, **de_random_kwargs(seed))
    x=np.asarray(result.x,dtype=float)
    d=diagnostics(x,target)
    row=dict(target_index=target_idx, restart=repeat,seed=seed,seed_targets=options['target_seed'],loss=float(result.fun),
             optimizer_success=bool(result.success),status='completed',nfev=int(result.nfev),nit=int(result.nit),
             mass_pass=d['mass_pass'],aux_pass=d['aux_pass'],all_pass=d['all_pass'],
             optimizer_constraints_ok=d['optimizer_constraints_ok'],
             target_H1=target[0],target_H2=target[1],target_H3=target[2],
             mass_H1=d['masses'][0],mass_H2=d['masses'][1],mass_H3=d['masses'][2],
             mc=d['mc'],kappa2V=d['kappa'], z_H1=d['z'][0],z_H2=d['z'][1],z_H3=d['z'][2],
             z_charged=d['zmc'],z_kappa=d['zk'],minimum_eigenvalue=d['evalmin'],message=str(result.message))
    row.update(dict(zip(PARAMETERS,x)))
    return row

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--preset',choices=['pilot','production'],default='pilot')
    p.add_argument('--targets',type=int,default=None)
    p.add_argument('--restarts',type=int,default=3)
    p.add_argument('--maxiter',type=int,default=None)
    p.add_argument('--popsize',type=int,default=None)
    p.add_argument('--target-seed',type=int,default=20261020)
    p.add_argument('--restart-seeds',type=int,nargs='*',default=None)
    p.add_argument('--output',type=Path,default=None)
    p.add_argument('--no-polish',action='store_true')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--limit-runs',type=int,default=None,help='Allow a bounded chunk for Colab sessions')
    a=p.parse_args()
    n=a.targets or (30 if a.preset=='pilot' else 1000)
    maxiter=a.maxiter or (150 if a.preset=='pilot' else 300)
    popsize=a.popsize or (12 if a.preset=='pilot' else 50)
    if n<1 or a.restarts<1:raise ValueError('positive targets/restarts required')
    seeds=a.restart_seeds or [11,23,37,47,59,71,83][:a.restarts]
    if len(seeds)!=a.restarts:raise ValueError('provide one seed per restart')
    out=a.output or Path(__file__).resolve().parents[1]/'data'/'new_runs'/a.preset
    out.mkdir(parents=True,exist_ok=True)
    output=out/'runs.csv'; targetfile=out/'targets.csv'; configfile=out/'config.json'
    config=dict(preset=a.preset,n_targets=n,n_restarts=a.restarts,restart_seeds=seeds,target_seed=a.target_seed,
                maxiter=maxiter,popsize=popsize,tol=.01,polish=not a.no_polish,v_ew=V_EW,
                mass_centers=MASS_CENTERS.tolist(),mass_widths=MASS_WIDTHS.tolist(),
                optimizer_charged_mass_squared_residual=[-260,260],optimizer_kappa_range=[.98,1.08],
                analysis_z_threshold=1,python=platform.python_version(),scipy=scipy.__version__,numpy=np.__version__,
                version_note='Fresh numerical experiment; never merge with legacy results.')
    if configfile.exists():
        stored=json.loads(configfile.read_text())
        for k in ['n_targets','n_restarts','restart_seeds','target_seed','maxiter','popsize','polish','analysis_z_threshold']:
            if stored[k]!=config[k]:raise ValueError(f'Existing config conflict on {k}: {stored[k]} vs {config[k]}. Choose new --output.')
    else:configfile.write_text(json.dumps(config,indent=2))
    rng=np.random.default_rng(a.target_seed)
    targets=rng.normal(MASS_CENTERS,MASS_WIDTHS,size=(n,3))
    if targetfile.exists():
        loaded=np.loadtxt(targetfile,delimiter=',',skiprows=1,ndmin=2)
        np.testing.assert_array_equal(loaded,targets)
    else:
        np.savetxt(targetfile,targets,delimiter=',',header='target_H1,target_H2,target_H3',comments='',fmt='%.17g')
    existing=[]
    if output.exists():
        if not a.resume:raise ValueError(f'{output} already exists. Use --resume')
        with output.open() as f:existing=list(csv.DictReader(f))
    completed={(int(r['target_index']),int(r['restart'])) for r in existing}
    done_now=0
    try:
        for idx,target in enumerate(targets):
            for repeat,base_seed in enumerate(seeds):
                if (idx,repeat) in completed:continue
                # Per-target seeds ensure that interrupted/resumed execution is deterministic.
                seed=int(base_seed+100003*idx)
                row=run_one(target,idx,repeat,seed,config)
                existing.append(row)
                atomic_write_rows(output,existing)
                done_now+=1
                print(f'target={idx+1}/{n} restart={repeat+1}/{a.restarts}: loss={row["loss"]:.4g}, '
                      f'zmax={max(abs(float(row[f"z_H{i}"])) for i in range(1,4)):.2f}, '
                      f'accepted={row["all_pass"]}, nfev={row["nfev"]}',flush=True)
                if a.limit_runs and done_now>=a.limit_runs:
                    print('Stopped after bounded chunk; resume with --resume');return
    except KeyboardInterrupt:
        print('Interrupted; completed rows saved. Use --resume',file=sys.stderr)
        raise
    print(f'All {n*a.restarts} runs complete. Dataset {output}')

if __name__=='__main__':main()
