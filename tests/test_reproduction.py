from pathlib import Path
import sys
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from model import mass_matrix,spectrum,V_EW,PARAMETERS

def test_reproduction_and_selection():
    df=pd.read_csv(ROOT/'data/raw/runs.csv')
    assert len(df)==3000
    assert df[['target_index','restart']].drop_duplicates().shape[0]==3000
    z=np.abs(df[['z_H1','z_H2','z_H3','z_charged','z_kappa']].to_numpy(float))
    assert (z<1).all(axis=1).sum()==2971
    assert df.loc[(z<1).all(axis=1),'target_index'].nunique()==1000

def test_real_numerical_example():
    df=pd.read_csv(ROOT/'data/raw/runs.csv')
    row=df.loc[(df.target_index==0)&(df.restart==2)].iloc[0]
    x=row[list(PARAMETERS)].to_numpy(float)
    M=mass_matrix(x)
    np.testing.assert_allclose(M,M.T,atol=1e-10)
    m,R,k,charged,_=spectrum(x)
    np.testing.assert_allclose(R.T@M@R,np.diag(m*m),atol=1e-8)
    np.testing.assert_allclose(m,row[['mass_H1','mass_H2','mass_H3']].to_numpy(float),atol=1e-8)
    assert np.max(np.abs((m-row[['target_H1','target_H2','target_H3']].to_numpy(float))/np.array([.05,.09,.09])))<1
    np.testing.assert_allclose([charged,k],[row.mc,row.kappa2V],atol=1e-8)

def test_fast_kernel_matches_reference():
    # If Numba is available, this also exercises compilation; otherwise it tests
    # the pure-Python fallback used by the optimizer.
    from fast_model import objective_fast, constraint_fast
    from model import spectral_loss, constraints
    df=pd.read_csv(ROOT/'data/raw/runs.csv').head(3)
    for _,row in df.iterrows():
        x=row[list(PARAMETERS)].to_numpy(float)
        target=row[['target_H1','target_H2','target_H3']].to_numpy(float)
        np.testing.assert_allclose(objective_fast(x,target),spectral_loss(x,target),rtol=1e-5,atol=1e-9)
        np.testing.assert_allclose(constraint_fast(x),constraints(x),rtol=1e-9,atol=1e-8)
