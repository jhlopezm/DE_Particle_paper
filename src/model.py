"""2HDM+S cubic CP-even forward model, with column-eigenvector convention.

R.T @ M2 @ R = diag(mass**2), phi = R @ H, H = R.T @ phi.
Important: sign choice for each eigenstate is arbitrary; H2 is oriented such
that kappa_2V >= 0.  All masses and VEVs in GeV.
"""
import numpy as np

V_EW = 246.22
PARAMETERS = ("lambda_1", "lambda_2", "lambda_3", "lambda_4", "lambda_sigma", "lambda_1sigma", "lambda_2sigma", "v_1", "v_sigma", "mu")
BOUNDS = ((0,1.5),(0,1.5),(-1.5,1.5),(-1.5,1.5),(0,1.5),(-1.5,1.5),(-1.5,1.5),(.01,4.),(250.,2000.),(-77.3,0.))
MASS_CENTERS = np.array([95.40, 125.09, 151.50])
MASS_WIDTHS = np.array([.05,.09,.09])
MC_TARGET, MC_WIDTH = 130., 1.
KAPPA_TARGET, KAPPA_WIDTH = 1.03, .05

def mass_matrix(x, v=V_EW):
    l1,l2,l3,l4,ls,l1s,l2s,v1,vs,mu = np.asarray(x,dtype=float)
    if not (0 < v1 < v and vs > 0):
        raise ValueError('Unphysical vacuum expectation values')
    v2=np.sqrt(v*v-v1*v1)
    m11=2*l1*v1*v1-mu*v2*vs/(np.sqrt(2)*v1)
    m12=mu*vs/np.sqrt(2)+v1*v2*(l3+l4)
    m13=mu*v2/np.sqrt(2)+l1s*v1*vs
    m22=2*l2*v2*v2-mu*v1*vs/(np.sqrt(2)*v2)
    m23=mu*v1/np.sqrt(2)+l2s*v2*vs
    m33=2*ls*vs*vs-mu*v1*v2/(np.sqrt(2)*vs)
    return np.array([[m11,m12,m13],[m12,m22,m23],[m13,m23,m33]],dtype=float)

def spectrum(x):
    """Return (m1,m2,m3), R, kappav, Mc, min eigenvalue.
    For a negative eigenvalue the corresponding mass is NaN (not zero).
    """
    evals,R=np.linalg.eigh(mass_matrix(x))
    masses = np.where(evals >= 0, np.sqrt(np.maximum(evals,0)), np.nan)
    v1=float(x[7]); v2=np.sqrt(V_EW**2-v1*v1)
    overlap=(v1*R[0,1]+v2*R[1,1])/V_EW
    if overlap < 0: R[:,1]*=-1
    kappa=abs(overlap)
    charged2=-V_EW**2/2*(np.sqrt(2)*x[9]*x[8]/(v1*v2)+x[3])
    charged=np.sqrt(charged2) if charged2 >= 0 else np.nan
    return masses,R,kappa,charged,evals[0]

def constraints(x):
    """Optimizer constraint residuals: (Mc^2-130^2) in GeV^2, kappa_2V.
    The original run constrained the squared-mass residual to +/-260 GeV^2.
    """
    v1=x[7]; v2=np.sqrt(V_EW**2-v1**2)
    mc2=-V_EW**2/2*(np.sqrt(2)*x[9]*x[8]/(v1*v2)+x[3])
    eig,R=np.linalg.eigh(mass_matrix(x))
    kappa=abs((v1*R[0,1]+v2*R[1,1])/V_EW)
    return np.array([mc2-MC_TARGET**2,kappa])

def spectral_loss(x,target):
    evals=np.linalg.eigvalsh(mass_matrix(x))
    if np.any(evals <= 0):
        return 1e8 + 1e-3*float(np.sum(np.minimum(evals,0)**2))
    masses=np.sqrt(evals)
    residual=(masses-np.asarray(target))/MASS_WIDTHS
    return float(residual @ residual)

def diagnostics(x,target):
    masses,R,kappa,mc,emin=spectrum(x)
    z=(masses-np.asarray(target))/MASS_WIDTHS
    zmc=abs(mc-MC_TARGET)/MC_WIDTH
    zk=abs(kappa-KAPPA_TARGET)/KAPPA_WIDTH
    goodmass=bool(np.all(np.isfinite(z)) and np.all(abs(z)<1))
    # Final uniform one-width post-run acceptance, distinct from optimizer bounds.
    goodaux=bool(np.isfinite(zmc) and zmc<1 and zk<1)
    strict=bool(np.isfinite(mc) and abs(mc*mc-130**2) <= 260 and .98 <= kappa <= 1.08)
    return dict(masses=masses,kappa=kappa,mc=mc,evalmin=emin,z=z,zmc=zmc,zk=zk,
                mass_pass=goodmass,aux_pass=goodaux,all_pass=goodmass and goodaux,
                optimizer_constraints_ok=strict)
