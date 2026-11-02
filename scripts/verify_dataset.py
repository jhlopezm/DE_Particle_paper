"""Verify immutable production data and final five-one-width acceptance.

Run from any working directory:
    python scripts/verify_dataset.py
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from model import PARAMETERS, MASS_WIDTHS, spectrum

r=pd.read_csv(ROOT/'data/raw/runs.csv')
t=pd.read_csv(ROOT/'data/raw/targets.csv')
assert len(r)==3000 and len(t)==1000
assert not r[['target_index','restart']].duplicated().any()
assert r.target_index.nunique()==1000
assert (r.groupby('target_index').size()==3).all()
np.testing.assert_allclose(r[['target_H1','target_H2','target_H3']].to_numpy(),t.to_numpy()[r.target_index.to_numpy()],rtol=0,atol=1e-11)

par=r[list(PARAMETERS)].to_numpy(float)
y=r[['target_H1','target_H2','target_H3']].to_numpy(float)
forward=np.array([np.r_[spectrum(x)[0],spectrum(x)[3],spectrum(x)[2]] for x in par])
stored=r[['mass_H1','mass_H2','mass_H3','mc','kappa2V']].to_numpy(float)
np.testing.assert_allclose(forward,stored,rtol=0,atol=1e-8)
z=np.column_stack([np.abs((forward[:,:3]-y)/MASS_WIDTHS),np.abs((forward[:,3]-130)/1),np.abs((forward[:,4]-1.03)/.05)])
orig=np.abs(r[['z_H1','z_H2','z_H3','z_charged','z_kappa']].to_numpy(float))
np.testing.assert_allclose(z,orig,rtol=0,atol=1e-7)
mask=(z<1).all(axis=1)
assert np.all(np.isfinite(z)) and mask.sum()==2971
assert ((z[:,:3]<1).all(axis=1)).sum()==2977
assert (z[:,3]<1).sum()==2994
assert (z[:,4]<1).sum()==3000
assert r.loc[mask,'target_index'].nunique()==1000
print('Verified 3000 unique optimizations, 1000 targets and 3 restarts per target.')
print('Final 5 x 1-sigma selection:',mask.sum(),'/ 3000; covered targets:',r.loc[mask,'target_index'].nunique())
print('Max absolute difference in reconstructed observables:',np.max(np.abs(forward-stored)))
