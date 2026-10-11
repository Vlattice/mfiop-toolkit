import numpy as np, re
HESS="/root/.claude/uploads/4547bb7a-db79-55d1-a64f-b9f7ee1ff460/b0d8c0b1-skatole_optfreq.hess"
txt=open(HESS).read().splitlines()
def find(sec):
    for i,l in enumerate(txt):
        if l.strip()==sec: return i
    return -1
# atoms
ia=find("$atoms"); nat=int(txt[ia+1])
syms=[]; mass=[]; coord=[]
for k in range(nat):
    p=txt[ia+2+k].split(); syms.append(p[0]); mass.append(float(p[1]))
    coord.append([float(p[2]),float(p[3]),float(p[4])])  # Bohr
coord=np.array(coord); mass=np.array(mass)
N=nat*3
# normal modes matrix
im=find("$normal_modes")
dims=txt[im+1].split(); nr,nc=int(dims[0]),int(dims[1])
M=np.zeros((nr,nc)); r=im+2
# ORCA prints in column blocks of 5
col=0
while col<nc:
    header=txt[r].split(); r+=1
    cols=[int(x) for x in header]
    for row in range(nr):
        parts=txt[r].split(); r+=1
        ridx=int(parts[0])
        for j,c in enumerate(cols):
            M[ridx,c]=float(parts[1+j])
    col=cols[-1]+1
# vibrational freqs
iv=find("$vibrational_frequencies"); nvf=int(txt[iv+1]); freqs=[]
for k in range(nvf):
    freqs.append(float(txt[iv+2+k].split()[1]))
freqs=np.array(freqs)
# molecular plane via ring heavy atoms (first ~9 C + N excluding methyl C which is index0)
# heavy atoms indices (non-H)
heavy=[i for i,s in enumerate(syms) if s!='H']
Hv=coord[heavy]; c0=Hv.mean(0); u,s,vt=np.linalg.svd(Hv-c0); nrm=vt[2]  # plane normal
print("plane normal:", np.round(nrm,3))
print(f"{'mode':>4} {'freq':>9} {'out-of-plane%':>13} {'topC_atoms(disp)':>18}")
for m in [25,26,27,28]:
    v=M[:,m].reshape(nat,3)   # cartesian displacement per atom
    # out of plane fraction: component along plane normal
    oop=np.array([np.dot(v[i],nrm) for i in range(nat)])
    ip = v - np.outer(oop,nrm)
    oop_frac=np.sum(oop**2)/np.sum(v**2)
    amp=np.linalg.norm(v,axis=1)
    top=np.argsort(-amp)[:4]
    print(f"{m:>4} {freqs[m]:9.2f} {100*oop_frac:12.1f}%  "+ " ".join(f"{syms[t]}{t}:{amp[t]:.2f}" for t in top))
np.save("/tmp/rework/nm.npy", M); 
np.save("/tmp/rework/coord_bohr.npy", coord)
import json; json.dump({"syms":syms,"mass":mass.tolist(),"freqs":freqs.tolist()}, open("/tmp/rework/meta.json","w"))
print("\nsaved nm.npy, coord_bohr.npy, meta.json")
