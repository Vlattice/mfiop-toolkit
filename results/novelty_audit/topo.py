"""
Brazo 2 (control topologico): red de proximidad de los Trp.
Pregunta: la separacion 5SYF vs 3J6F que da MFIOP, la da tambien la pura topologia?
Metricas: espectro de la Laplaciana, grado medio, clustering, Fiedler value.
Y el test clave: la autoestructura de MFIOP (con orientaciones) aporta algo que
la matriz de adyacencia PURA (solo distancias, SIN orientaciones) no da?
"""
import numpy as np, networkx as nx
from Bio.PDB import PDBParser
from itertools import combinations

def load_trp(path):
    p=PDBParser(QUIET=True); model=list(p.get_structure("x",path).get_models())[0]
    pos=[]; dip=[]
    for ch in model:
        for res in ch:
            if res.get_resname()!="TRP": continue
            try:
                cg=res["CG"].get_vector().get_array(); ce2=res["CE2"].get_vector().get_array()
            except KeyError: continue
            ax=ce2-cg; n=np.linalg.norm(ax)
            if n<1e-6: continue
            pos.append(cg); dip.append(ax/n)
    return np.array(pos,float), np.array(dip,float)

def topo_metrics(pos, cutoff=25.0):
    N=len(pos); A=np.zeros((N,N))
    for i,j in combinations(range(N),2):
        d=np.linalg.norm(pos[j]-pos[i])
        if d<=cutoff: A[i,j]=A[j,i]=1.0
    G=nx.from_numpy_array(A)
    deg=A.sum(1)
    # Laplaciana y su espectro
    L=np.diag(deg)-A
    evL=np.sort(np.linalg.eigvalsh(L))
    fiedler=evL[1] if N>1 else 0.0   # algebraic connectivity
    # adjacency ponderada por 1/r^3 SIN orientacion (topologia geometrica pura)
    W=np.zeros((N,N))
    for i,j in combinations(range(N),2):
        d=np.linalg.norm(pos[j]-pos[i])
        if d>1e-6: W[i,j]=W[j,i]=1.0/d**3
    evW=np.linalg.eigvalsh(W)
    return dict(N=N, meandeg=float(deg.mean()),
                clustering=float(nx.average_clustering(G)),
                fiedler=float(fiedler),
                W_range=float(evW.max()-evW.min()))

def mfiop_range(pos,dip):
    N=len(pos); J=np.zeros((N,N))
    for i,j in combinations(range(N),2):
        r=pos[j]-pos[i]; d=np.linalg.norm(r)
        if d<1e-6: continue
        rh=r/d; kap=dip[i]@dip[j]-3*(dip[i]@rh)*(dip[j]@rh)
        J[i,j]=J[j,i]=kap/d**3
    ev=np.linalg.eigvalsh(J)
    return float(ev.max()-ev.min())

print(f"{'struct':>6} {'N':>4} {'meandeg':>8} {'clust':>6} {'fiedler':>8} {'W_range(noOri)':>14} {'MFIOP_range(Ori)':>16}")
rows={}
for name,path in [("5SYF","5SYF.pdb"),("3J6F","3J6F.pdb")]:
    pos,dip=load_trp(path); t=topo_metrics(pos); m=mfiop_range(pos,dip)
    rows[name]=(t,m)
    print(f"{name:>6} {t['N']:>4} {t['meandeg']:>8.2f} {t['clustering']:>6.2f} {t['fiedler']:>8.3f} {t['W_range']:>14.3e} {m:>16.3e}")

print("\n--- TEST CLAVE: MFIOP aporta algo mas que la geometria sin orientacion? ---")
print("Comparo el ratio 3J6F/5SYF en cada descriptor. Si MFIOP (con orientacion)")
print("da un ratio DISTINTO al de W_range (sin orientacion), la orientacion aporta.")
t5,m5=rows["5SYF"]; t3,m3=rows["3J6F"]
print(f"  ratio W_range (sin orient) 3J6F/5SYF = {t3['W_range']/t5['W_range']:.3f}")
print(f"  ratio MFIOP  (con orient)  3J6F/5SYF = {m3/m5:.3f}")
print(f"  ratio Fiedler (topol pura) 3J6F/5SYF = {t3['fiedler']/t5['fiedler']:.3f}" if t5['fiedler']>0 else "  fiedler 5SYF=0")
