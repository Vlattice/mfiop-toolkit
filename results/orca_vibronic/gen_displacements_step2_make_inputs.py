import numpy as np, json, os
M=np.load("nm.npy"); coordB=np.load("coord_bohr.npy"); meta=json.load(open("meta.json"))
syms=meta["syms"]; mass=np.array(meta["mass"]); freqs=np.array(meta["freqs"])
nat=len(syms); BOHR=0.52917721067
X0=coordB*BOHR   # optimized geometry in Angstrom
OUT="/mnt/user-data/outputs/ORCA_Carlos_ronda2"

def xyzblock(X):
    return "\n".join(f"{syms[i]:2s} {X[i,0]:14.8f} {X[i,1]:14.8f} {X[i,2]:14.8f}" for i in range(nat))

HEADER="""! CAM-B3LYP def2-SVP D3BJ TightSCF
%maxcore 3000
%pal nprocs 1 end
%tddft
  nroots 6
  tda false
end
* xyz 0 1
"""
def write_inp(name, X, comment):
    with open(f"{OUT}/{name}.inp","w") as f:
        f.write(f"# {comment}\n# correr EN SERIE (nprocs 1). TD-DFT para el dipolo de transicion 1La.\n")
        f.write(HEADER); f.write(xyzblock(X)+"\n*\n")

manifest={"BOHR":BOHR,"modes":{}}
# equilibrium
write_inp("skatole_td_eq", X0, "Geometria de equilibrio (optimizada). 1La de referencia.")
STEP_ANG=0.10
for m in [26,28]:
    L=M[:,m].reshape(nat,3)
    norm_w=np.sqrt(np.sum(mass[:,None]*L**2))
    disp_dir=L/norm_w                       # Cartesian displacement per unit Q (Ang per amu^0.5-Ang)
    per_atom=np.linalg.norm(disp_dir,axis=1)
    dQ=STEP_ANG/per_atom.max()              # amu^0.5 * Ang, so that max atomic step = 0.10 Ang
    dX=disp_dir*dQ
    Xp=X0+dX; Xm=X0-dX
    write_inp(f"skatole_td_m{m}_plus",  Xp, f"Modo {m} ({freqs[m]:.1f} cm-1), desplazado +dQ")
    write_inp(f"skatole_td_m{m}_minus", Xm, f"Modo {m} ({freqs[m]:.1f} cm-1), desplazado -dQ")
    manifest["modes"][str(m)]={"freq_cm":float(freqs[m]),"norm_w":float(norm_w),
        "dQ_amu05_Ang":float(dQ),"max_atom_step_Ang":float(per_atom.max()*dQ)}
    print(f"mode {m}: freq={freqs[m]:.2f}  dQ={dQ:.5f} amu^0.5*Ang  maxstep={per_atom.max()*dQ:.3f} Ang")
json.dump(manifest, open(f"{OUT}/manifest_round2.json","w"), indent=2)
print("\nfiles:"); print("\n".join(sorted(os.listdir(OUT))))
