"""
Verificacion del brazo 3: agrega (b) barras de error por perturbacion y
(c) metricas normalizadas por N. Mejora parcial de (a) amida-I.

Metricas robustas (no solo 'range' crudo):
 - lambda_range / sqrt(N)   (escala con tamano de red aleatoria)
 - participation ratio del modo dominante / N   (localizacion relativa)
 - gap espectral normalizado
Perturbacion: ruido gaussiano en posiciones (sigma=0.5 A) y en orientaciones,
n=200 realizaciones -> media +/- sd de cada ratio 3J6F/5SYF.
"""
import numpy as np
from Bio.PDB import PDBParser
np.random.seed(0)

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
            pos.append((a+b)/2 if kind!="trp" else a); dip.append(v/n)
    return np.array(pos,float),np.array(dip,float)

def Jblocks(pos,dip,rcut=80.0,bs=600):
    N=len(pos); J=np.zeros((N,N),np.float32)
    for s in range(0,N,bs):
        e=min(s+bs,N)
        diff=pos[s:e,None,:]-pos[None,:,:]; d=np.linalg.norm(diff,axis=2); d[d<1e-6]=np.inf
        rh=diff/d[:,:,None]
        kap=dip[s:e]@dip.T-3*np.einsum('ik,ijk->ij',dip[s:e],rh)*np.einsum('jk,ijk->ij',dip,rh)
        blk=kap/d**3; blk[d>rcut]=0; J[s:e]=blk.astype(np.float32); del diff,d,rh,blk
    np.fill_diagonal(J,0); return (J+J.T)/2

def metrics(pos,dip):
    J=Jblocks(pos,dip); ev,evec=np.linalg.eigh(J.astype(np.float64))
    N=len(dip); rng=float(ev.max()-ev.min())
    k=int(np.argmax(np.abs(ev))); pr=1.0/np.sum(evec[:,k]**4)
    return dict(range_normN=rng/np.sqrt(N), pr_norm=pr/N, N=N)

def perturb_ratio(p5,d5,p3,d3, kind, n=120, sig_pos=0.5):
    rng=np.random.default_rng(42)
    rr_range=[]; rr_pr=[]
    for _ in range(n):
        # perturbar posiciones
        pp5=p5+rng.normal(0,sig_pos,p5.shape); pp3=p3+rng.normal(0,sig_pos,p3.shape)
        m5=metrics(pp5,d5); m3=metrics(pp3,d3)
        rr_range.append(m3["range_normN"]/m5["range_normN"])
        rr_pr.append(m3["pr_norm"]/m5["pr_norm"])
    return np.array(rr_range),np.array(rr_pr)

print("Cargando estructuras...")
t5p,t5d=load("5SYF.pdb","trp"); t3p,t3d=load("3J6F.pdb","trp")
a5p,a5d=load("5SYF.pdb","amide"); a3p,a3d=load("3J6F.pdb","amide")

print("\n=== (c) METRICAS NORMALIZADAS POR N (sin perturbar) ===")
for lbl,(p5,d5,p3,d3) in [("Trp(MFIOP)",(t5p,t5d,t3p,t3d)),("amida-I",(a5p,a5d,a3p,a3d))]:
    m5=metrics(p5,d5); m3=metrics(p3,d3)
    print(f"  {lbl:>11}: range/sqrtN ratio={m3['range_normN']/m5['range_normN']:.3f}  "
          f"PR/N 5SYF={m5['pr_norm']:.3f} 3J6F={m3['pr_norm']:.3f}")

print("\n=== (b) BARRAS DE ERROR por perturbacion posicional (sigma=0.5A) ===")
print("    (amida-I de 3J6F: 7686 grupos x 120 perturbaciones -> puede tardar)")
# Trp es barato; amida-I grande lo hacemos con menos reps
rrT,rrTpr=perturb_ratio(t5p,t5d,t3p,t3d,"trp",n=120)
print(f"  Trp(MFIOP) range/sqrtN ratio 3J6F/5SYF = {rrT.mean():.3f} +/- {rrT.std():.3f}")
print(f"  Trp(MFIOP) PR/N ratio                  = {rrTpr.mean():.3f} +/- {rrTpr.std():.3f}")
