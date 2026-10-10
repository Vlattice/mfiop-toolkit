"""Brazo 3 - amida-I por bloques (memoria controlada) para estructuras grandes."""
import numpy as np
from Bio.PDB import PDBParser
np.random.seed(0)

def load_amide(path):
    p=PDBParser(QUIET=True); model=list(p.get_structure("x",path).get_models())[0]
    pos=[]; dip=[]
    for ch in model:
        for res in ch:
            try:
                C=res["C"].get_vector().get_array(); O=res["O"].get_vector().get_array()
            except KeyError: continue
            co=O-C; n=np.linalg.norm(co)
            if n<1e-6: continue
            pos.append((C+O)/2.0); dip.append(co/n)
    return np.array(pos,float), np.array(dip,float)

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

def build_J_blocks(pos,dip,rcut=80.0,bs=600):
    N=len(pos); J=np.zeros((N,N),dtype=np.float32)
    for a in range(0,N,bs):
        b=min(a+bs,N)
        diff=pos[a:b,None,:]-pos[None,:,:]       # (bs,N,3)
        d=np.linalg.norm(diff,axis=2)            # (bs,N)
        d[d<1e-6]=np.inf
        rh=diff/d[:,:,None]
        eiej=dip[a:b]@dip.T
        eir=np.einsum('ik,ijk->ij',dip[a:b],rh)
        ejr=np.einsum('jk,ijk->ij',dip,rh.transpose(1,0,2)).T if False else np.einsum('jk,ijk->ij',dip,rh)
        kap=eiej-3*eir*ejr
        blk=kap/(d**3); blk[d>rcut]=0.0
        J[a:b,:]=blk.astype(np.float32)
        del diff,d,rh,blk
    np.fill_diagonal(J,0.0); J=(J+J.T)/2
    return J

def obs(J):
    ev=np.linalg.eigvalsh(J.astype(np.float64))
    return float(ev.max()-ev.min())

def obs_full(J):
    ev,evec=np.linalg.eigh(J.astype(np.float64))
    k=int(np.argmax(np.abs(ev))); pr=1.0/np.sum(evec[:,k]**4)
    return float(ev.max()-ev.min()), float(pr)

R={}
for name,path in [("5SYF","5SYF.pdb"),("3J6F","3J6F.pdb")]:
    ap,ad=load_amide(path)
    Ja=build_J_blocks(ap,ad)
    ar=obs(Ja); del Ja
    tp,td=load_trp(path); Jt=build_J_blocks(tp,td); tr,tpr=obs_full(Jt)
    R[name]=(len(ap),ar,len(tp),tr,tpr)
    print(f"{name}: amideI N={len(ap)} range={ar:.3e} | Trp N={len(tp)} range={tr:.3e} PR={tpr:.2f}")

(a5N,a5r,t5N,t5r,t5pr)=R["5SYF"]; (a3N,a3r,t3N,t3r,t3pr)=R["3J6F"]
print("\n--- DECISIVO ---")
print(f"  amida-I range ratio 3J6F/5SYF    = {a3r/a5r:.3f}")
print(f"  Trp(MFIOP) range ratio 3J6F/5SYF = {t3r/t5r:.3f}")
print(f"  -> amida-I ya separa las estructuras con ese ratio.")
print(f"     Si el ratio Trp es similar, MFIOP es REDUNDANTE con amida-I.")
print(f"     Si es claramente distinto, Trp aporta eje independiente.")
