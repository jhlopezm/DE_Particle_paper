"""Optional Numba-accelerated optimizer kernel, equivalent to model.py.
For reproducibility the analysis always recomputes observables with model.py.
"""
import numpy as np
try:
    from numba import njit
    NUMBA_AVAILABLE = True
except ImportError:
    # When Numba is not installed, preserve the @njit(cache=True)
    # decorator syntax while executing the functions as ordinary Python.
    # This mode is correct but noticeably slower for large DE scans.
    NUMBA_AVAILABLE = False
    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]):
            return args[0]  # supports bare @njit
        def decorate(fn):
            return fn      # supports @njit(cache=True), etc.
        return decorate
from model import V_EW, MASS_WIDTHS, MASS_CENTERS

@njit(cache=True)
def _matrix(x):
    l1,l2,l3,l4,ls,l1s,l2s,v1,vs,mu=x
    v2=np.sqrt(V_EW*V_EW-v1*v1)
    m11=2*l1*v1*v1-mu*v2*vs/(np.sqrt(2.)*v1)
    m12=mu*vs/np.sqrt(2.)+v1*v2*(l3+l4)
    m13=mu*v2/np.sqrt(2.)+l1s*v1*vs
    m22=2*l2*v2*v2-mu*v1*vs/(np.sqrt(2.)*v2)
    m23=mu*v1/np.sqrt(2.)+l2s*v2*vs
    m33=2*ls*vs*vs-mu*v1*v2/(np.sqrt(2.)*vs)
    return np.array([[m11,m12,m13],[m12,m22,m23],[m13,m23,m33]])

@njit(cache=True)
def objective_fast(x,target):
    evals=np.linalg.eigvalsh(_matrix(x))
    if evals[0] <= 0:
        val=0.0
        for i in range(3):
            if evals[i] < 0: val += evals[i]*evals[i]
        return 1e8+1e-3*val
    loss=0.0
    for i in range(3):
        z=(np.sqrt(evals[i])-target[i])/MASS_WIDTHS[i]
        loss+=z*z
    return loss

@njit(cache=True)
def constraint_fast(x):
    v1,vs,mu=x[7],x[8],x[9]
    v2=np.sqrt(V_EW*V_EW-v1*v1)
    mc2=-V_EW*V_EW/2*(np.sqrt(2.)*mu*vs/(v1*v2)+x[3])
    evals,R=np.linalg.eigh(_matrix(x))
    kapp=abs((v1*R[0,1]+v2*R[1,1])/V_EW)
    return np.array([mc2-130.**2,kapp])
