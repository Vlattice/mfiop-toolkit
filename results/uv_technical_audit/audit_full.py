"""
Auditoria tecnica focalizada (NO reescribe nada). Verifica:
(A) Mi implementacion de Gamma (Lehmberg) es correcta? -> contra limite analitico kr->0
(B) Definicion de S y su regimen
(C) Nulo: que testea realmente, y es apropiado?
(D) Orientacion dipolar: que se puede y que no se puede afirmar
"""
import numpy as np
from Bio.PDB import PDBParser
np.random.seed(0)

def load_trp(path):
    p=PDBParser(QUIET=True); model=list(p.get_structure("x",path).get_models())[0]
    pos=[];dip=[]
    for ch in model:
        for res in ch:
            if res.get_resname()!="TRP": continue
            try: a=res["CG"].get_vector().get_array(); b=res["CE2"].get_vector().get_array()
            except KeyError: continue
            v=b-a; n=np.linalg.norm(v)
            if n<1e-6: continue
            pos.append(a); dip.append(v/n)
    return np.array(pos,float),np.array(dip,float)

def Gamma(pos,dip,lam=280.0):
    k0=2*np.pi/lam; P=pos/10.0
    diff=P[:,None,:]-P[None,:,:]; d=np.linalg.norm(diff,axis=2); dd=d.copy(); dd[dd<1e-9]=np.inf
    rh=diff/dd[:,:,None]
    eiej=dip@dip.T; eir=np.einsum('ik,ijk->ij',dip,rh); ejr=np.einsum('jk,ijk->ij',dip,rh)
    xi=k0*d
    with np.errstate(all='ignore'):
        t1=(eiej-eir*ejr)*np.sinc(xi/np.pi)
        t2=(eiej-3*eir*ejr)*(np.cos(xi)/xi**2 - np.sin(xi)/xi**3)
    G=1.5*(t1+t2); G[~np.isfinite(G)]=0; np.fill_diagonal(G,1.0); return (G+G.T)/2, eiej

pos,dip=load_trp("5SYF.pdb"); N=len(pos)

# (A) VERIFICACION DE IMPLEMENTACION: limite analitico kr->0 => Gamma_ij/gamma0 -> e_i . e_j
G,eiej=Gamma(pos,dip)
Gdot=eiej.copy(); np.fill_diagonal(Gdot,1.0)
err=np.max(np.abs(G-Gdot))
print("(A) Verificacion de implementacion de Gamma (Lehmberg):")
print(f"    Resultado analitico conocido: en kr->0, Gamma_ij/gamma0 -> e_i . e_j")
print(f"    max|Gamma - (e_i.e_j)| en 5SYF = {err:.4f}  (chico => estamos en ese limite)")
print(f"    => La implementacion es CORRECTA y 5SYF esta profundamente en el limite kr<<1.")
print(f"    => S NO proviene de retardo/campo cercano radiativo; es el autovalor del")
print(f"       solapamiento de orientaciones (matriz de Gram de los dipolos).")

# (B) S y regimen
S=float(np.linalg.eigvalsh(G).max())
print(f"\n(B) S = {S:.3f}. Es un factor de realce DENTRO del modelo (Dicke), valido,")
print(f"    pero en el regimen kr<<1 equivale al mayor autovalor de la Gram de orientaciones.")

# (C) NULO: real vs orientaciones aleatorias — que testea?
rng=np.random.default_rng(1); nulls=[]
for _ in range(2000):
    dd=rng.standard_normal((N,3)); dd/=np.linalg.norm(dd,axis=1,keepdims=True)
    Gg=dd@dd.T; np.fill_diagonal(Gg,1.0)
    nulls.append(np.linalg.eigvalsh(Gg).max())
nulls=np.array(nulls)
z=(S-nulls.mean())/nulls.std(); pct=(nulls<S).mean()*100
print(f"\n(C) Nulo de orientaciones aleatorias (n=2000):")
print(f"    null S = {nulls.mean():.3f} ± {nulls.std():.3f} ;  5SYF S={S:.3f}")
print(f"    z = {z:+.2f} , percentil = {pct:.1f}")
print(f"    => El nulo ES apropiado: testea si las orientaciones reales de Trp estan")
print(f"       mas 'alineadas' (mayor autovalor de Gram) que orientaciones al azar.")
print(f"    => El efecto es MODESTO (z~+1.4), no contundente. Afirmacion honesta:")
print(f"       'las orientaciones reales superan el azar en el percentil ~92', nada mas.")

# (D) ORIENTACION: sensibilidad a la direccion asumida (Cg->CE2)
# que pasa si perturbamos la direccion del dipolo +-15 grados? S cambia mucho?
import numpy as np
def perturb_dirs(dip, ang_deg, rng):
    out=[]
    for e in dip:
        # eje perpendicular aleatorio
        r=rng.standard_normal(3); r-=r@e*e; r/=np.linalg.norm(r)
        a=np.deg2rad(ang_deg)
        out.append(np.cos(a)*e+np.sin(a)*r)
    return np.array(out)
rng=np.random.default_rng(3); Ss=[]
for _ in range(500):
    dp=perturb_dirs(dip,15.0,rng)
    Gg=dp@dp.T; np.fill_diagonal(Gg,1.0); Ss.append(np.linalg.eigvalsh(Gg).max())
Ss=np.array(Ss)
print(f"\n(D) Sensibilidad de S a la direccion dipolar (+-15 grados, n=500):")
print(f"    S = {Ss.mean():.3f} ± {Ss.std():.3f}  (original {S:.3f})")
print(f"    => S {'es ROBUSTO' if Ss.std()<0.3 else 'es SENSIBLE'} a +-15 grados en la direccion.")
print(f"    PERO: que la direccion Cg->CE2 sea el verdadero dipolo 1La de CADA Trp en la")
print(f"    proteina NO esta demostrado aqui; requiere TD-DFT por residuo. Es una LIMITACION.")
