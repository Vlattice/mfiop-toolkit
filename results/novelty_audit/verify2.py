"""
Verificacion rigurosa. Arregla 3 cosas que marco la auditoria:
 1. amida-I tambien con barras de error (antes era un punto sin error).
 2. MISMA definicion de metrica y PR para ambos modelos (una sola funcion).
 3. Controla el confound de tamano: subsamplea la red grande a la N de la chica,
    con muchos draws -> comparacion justa Y barras de error naturales.

Definiciones explicitas (para auditar):
 - J_ij = kappa_ij / r_ij^3 ,  kappa = e_i.e_j - 3(e_i.r)(e_j.r)  [identica ambos]
 - participation ratio PR = 1 / sum_i v_i^4 , v = autovector del modo de MAYOR |lambda|
 - range = lambda_max - lambda_min
 - subsample: elegир N_sub residuos al azar (sin reemplazo); N_sub = min(N_5SYF, N_3J6F)
   para cada tipo. Trp: N_sub=8.  amida-I: N_sub=852.
 - n_draws realizaciones -> media y sd de cada metrica POR ESTRUCTURA.
"""
import numpy as np
from Bio.PDB import PDBParser

def load(path, kind):
    p=PDBParser(QUIET=True); model=list(p.get_structure("x",path).get_models())[0]
    pos=[];dip=[]
    for ch in model:
        for res in ch:
            if kind=="trp":
                if res.get_resname()!="TRP": continue
                try: a=res["CG"].get_vector().get_array(); b=res["CE2"].get_vector().get_array()
                except KeyError: continue
            else:
                try: a=res["C"].get_vector().get_array(); b=res["O"].get_vector().get_array()
                except KeyError: continue
            v=b-a; n=np.linalg.norm(v)
            if n<1e-6: continue
            pos.append(a if kind=="trp" else (a+b)/2); dip.append(v/n)
    return np.array(pos,float),np.array(dip,float)

def Jmat(pos,dip,rcut=80.0):
    diff=pos[:,None,:]-pos[None,:,:]; d=np.linalg.norm(diff,axis=2); 
    dd=d.copy(); dd[dd<1e-6]=np.inf; rh=diff/dd[:,:,None]
    kap=dip@dip.T-3*np.einsum('ik,ijk->ij',dip,rh)*np.einsum('jk,ijk->ij',dip,rh)
    J=kap/(dd**3); J[d>rcut]=0; np.fill_diagonal(J,0); return (J+J.T)/2

def met(pos,dip):
    J=Jmat(pos,dip); ev,evec=np.linalg.eigh(J)
    k=int(np.argmax(np.abs(ev)))
    return (ev.max()-ev.min()), 1.0/np.sum(evec[:,k]**4)

def subsample_stats(pos,dip,Nsub,ndraws=200,seed=1):
    rng=np.random.default_rng(seed); N=len(pos)
    rngs=[]; prs=[]
    for _ in range(ndraws):
        idx=rng.choice(N,Nsub,replace=False) if N>Nsub else np.arange(N)
        r,pr=met(pos[idx],dip[idx]); rngs.append(r); prs.append(pr)
    return np.array(rngs),np.array(prs)

print("Comparacion JUSTA: subsampleando la red grande a la N de la chica\n")
for kind in ["trp","amide"]:
    p5,d5=load("5SYF.pdb",kind); p3,d3=load("3J6F.pdb",kind)
    Nsub=min(len(p5),len(p3))
    nd = 200 if kind=="trp" else 40   # amida-I mas caro
    r5,pr5=subsample_stats(p5,d5,Nsub,ndraws=nd,seed=1)
    r3,pr3=subsample_stats(p3,d3,Nsub,ndraws=nd,seed=2)
    print(f"=== {kind.upper()}  (N_sub={Nsub}, draws={nd}) ===")
    print(f"  range  5SYF={r5.mean():.3e}±{r5.std():.1e}   3J6F={r3.mean():.3e}±{r3.std():.1e}")
    # separacion correcta: diferencia de medias / sd combinada
    sep=(r3.mean()-r5.mean())/np.sqrt(r5.std()**2+r3.std()**2+1e-30)
    print(f"    separacion (range) = {sep:+.2f} sigma")
    print(f"  PR     5SYF={pr5.mean():.2f}±{pr5.std():.2f}   3J6F={pr3.mean():.2f}±{pr3.std():.2f}")
    sepp=(pr3.mean()-pr5.mean())/np.sqrt(pr5.std()**2+pr3.std()**2+1e-30)
    print(f"    separacion (PR)    = {sepp:+.2f} sigma\n")
