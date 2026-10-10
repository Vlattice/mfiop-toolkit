"""
Brazo 1 (MFIOP puro) + Brazo 4 (nulo geometrico) sobre 5SYF y 3J6F.
Descriptor: autoestructura de la matriz de acoplamiento dipolar de Trp.
Observables: lambda_range, gap espectral s_max, participation ratio del modo dominante,
             anisotropia. Comparados contra nulo de orientacion (500 draws).
"""
import numpy as np
from Bio.PDB import PDBParser
from itertools import combinations
np.random.seed(0)

def load_trp(path):
    p=PDBParser(QUIET=True); st=p.get_structure("x",path)
    model=list(st.get_models())[0]
    pos=[]; dip=[]
    for ch in model:
        for res in ch:
            if res.get_resname()!="TRP": continue
            try:
                cg=res["CG"].get_vector().get_array()
                ce2=res["CE2"].get_vector().get_array()
            except KeyError: continue
            ax=ce2-cg; n=np.linalg.norm(ax)
            if n<1e-6: continue
            pos.append(cg); dip.append(ax/n)
    return np.array(pos,float), np.array(dip,float)

def build_J(pos, dip, rcut=80.0):
    N=len(pos); J=np.zeros((N,N))
    for i,j in combinations(range(N),2):
        r=pos[j]-pos[i]; d=np.linalg.norm(r)
        if d>rcut or d<1e-6: continue
        rh=r/d
        kap=dip[i]@dip[j]-3*(dip[i]@rh)*(dip[j]@rh)
        J[i,j]=J[j,i]=kap/d**3
    return J

def observables(J, dip):
    ev,evec=np.linalg.eigh(J)
    lam_range=float(ev.max()-ev.min())
    # participation ratio del modo de mayor |lambda|
    k=int(np.argmax(np.abs(ev)))
    v=evec[:,k]; pr=1.0/np.sum(v**4)   # IPR^-1, 1..N
    # anisotropia: proyeccion del modo dominante sobre z vs x,y
    z=dip[:,2]; fz=(evec.T@z)**2
    fx=(evec.T@dip[:,0])**2; fy=(evec.T@dip[:,1])**2
    kz=int(np.argmax(fz))
    aniso=float(fz[kz]/(fx[kz]+fy[kz]+1e-12))
    # gap espectral maximo entre autovalores consecutivos
    s=np.sort(ev); smax=float(np.max(np.diff(s))) if len(s)>1 else 0.0
    return dict(N=len(dip), lam_range=lam_range, pr=float(pr), aniso=aniso, smax=smax)

def null_orientation(pos, nreps=500):
    """Brazo 4: mantiene posiciones, randomiza orientaciones de dipolos."""
    rng=np.random.default_rng(1); N=len(pos)
    out=[]
    for _ in range(nreps):
        d=rng.standard_normal((N,3)); d/=np.linalg.norm(d,axis=1,keepdims=True)
        J=build_J(pos,d)
        out.append(observables(J,d)["lam_range"])
    return np.array(out)

for name,path in [("5SYF","5SYF.pdb"),("3J6F","3J6F.pdb")]:
    pos,dip=load_trp(path)
    J=build_J(pos,dip)
    obs=observables(J,dip)
    null=null_orientation(pos, nreps=300)
    z=(obs["lam_range"]-null.mean())/null.std()
    pct=float((null<obs["lam_range"]).mean()*100)
    print(f"\n=== {name}  (N_TRP={obs['N']}) ===")
    print(f"  lambda_range (real) = {obs['lam_range']:.4e}")
    print(f"  null orient  mean±sd= {null.mean():.4e} ± {null.std():.4e}")
    print(f"  z-score             = {z:+.2f}   percentil = {pct:.1f}")
    print(f"  participation ratio = {obs['pr']:.2f} (de {obs['N']})")
    print(f"  anisotropia A_z/(Ax+Ay) = {obs['aniso']:.3f}")
    print(f"  s_max (gap espectral)   = {obs['smax']:.4e}")
