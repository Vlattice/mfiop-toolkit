"""
AUDITORIA TECNICA de los 3 puntos criticos del paper UV, sobre 5SYF.
(1) La matriz UV (radiancia) comparte REALMENTE la eigenstructure de la matriz IR (acoplamiento)?
(3) Que significa S=5.06? Esta dominado por el producto escalar de orientaciones (limite kr->0)
    o por el termino de campo cercano real de la interaccion dipolar?
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

pos,dip=load_trp("5SYF.pdb"); N=len(pos)
lam=280.0; k0=2*np.pi/lam  # nm^-1
# distancias en Angstrom -> nm
P=pos/10.0
diff=P[:,None,:]-P[None,:,:]; d=np.linalg.norm(diff,axis=2); dd=d.copy(); dd[dd<1e-9]=np.inf
rh=diff/dd[:,:,None]
eiej=dip@dip.T
eir=np.einsum('ik,ijk->ij',dip,rh); ejr=np.einsum('jk,ijk->ij',dip,rh)

# --- Matriz IR: acoplamiento dipolo-dipolo J = kappa/r^3 ---
kap=eiej-3*eir*ejr
J_IR=kap/dd**3; np.fill_diagonal(J_IR,0); J_IR=(J_IR+J_IR.T)/2

# --- Matriz UV: radiancia de Lehmberg Gamma_ij/gamma0 ---
xi=k0*d  # kr adimensional
with np.errstate(divide='ignore',invalid='ignore'):
    t1=(eiej-eir*ejr)*np.sinc(xi/np.pi)                      # sin(xi)/xi
    t2=(eiej-3*eir*ejr)*(np.cos(xi)/xi**2 - np.sin(xi)/xi**3)
G=1.5*(t1+t2); 
G[~np.isfinite(G)]=0; np.fill_diagonal(G,1.0); G=(G+G.T)/2

print(f"5SYF: N={N} TRP ;  k0*r rango = [{xi[d>0].min():.3f}, {xi[d>0].max():.3f}]")

# ---- (3) S y su origen ----
evG=np.linalg.eigvalsh(G); S=float(evG.max())
print(f"\n(3) S = max eig(Gamma/gamma0) = {S:.3f}")
# matriz de SOLO producto escalar de orientaciones (limite kr->0 del termino dominante)
Gdot=eiej.copy(); np.fill_diagonal(Gdot,1.0)
S_dot=float(np.linalg.eigvalsh(Gdot).max())
print(f"    S de SOLO (e_i . e_j)              = {S_dot:.3f}")
print(f"    -> si S ~ S_dot, S refleja orientaciones, no campo cercano real")
print(f"    diferencia relativa = {abs(S-S_dot)/S*100:.1f}%")

# ---- (1) misma eigenstructure IR vs UV? ----
evIR,VIR=np.linalg.eigh(J_IR)
evUV,VUV=np.linalg.eigh(G)
# ordenar por |autovalor| desc
iIR=np.argsort(-np.abs(evIR)); iUV=np.argsort(-np.abs(evUV))
VIR=VIR[:,iIR]; VUV=VUV[:,iUV]
# overlap del modo dominante IR con todos los UV
dom_overlap=np.abs(VIR[:,0]@VUV[:,0])
# matriz de solapamiento |<vIR_i|vUV_j>|
Ov=np.abs(VIR.T@VUV)
# cuanto del modo dominante IR cae en el dominante UV
print(f"\n(1) Eigenstructure IR vs UV (modo dominante):")
print(f"    |<v_dom^IR | v_dom^UV>| = {dom_overlap:.3f}  (1=identicos, 0=ortogonales)")
print(f"    max overlap del v_dom^IR con cualquier modo UV = {Ov[0].max():.3f} (modo UV #{int(np.argmax(Ov[0]))})")
# alineacion global: traza de Ov^2 / N (1 = misma base)
align=float(np.sum(Ov**2)/N)  # siempre 1 si ambas ortonormales... mejor: subespacio dominante
# alineacion de los 3 modos dominantes (subespacios)
k=3
PIR=VIR[:,:k]@VIR[:,:k].T; PUV=VUV[:,:k]@VUV[:,:k].T
subspace_align=float(np.trace(PIR@PUV)/k)  # 1 si los 3-dim subespacios coinciden
print(f"    alineacion subespacio 3 modos dominantes = {subspace_align:.3f} (1=mismo subespacio)")
