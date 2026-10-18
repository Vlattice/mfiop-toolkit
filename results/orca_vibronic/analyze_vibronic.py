import numpy as np, json, pickle
data=pickle.load(open("td.pkl","rb"))
man=json.load(open("/mnt/user-data/outputs/ORCA_Carlos_ronda2/manifest_round2.json"))
AU2D=2.5417464
# constants
hbar=1.054571817e-34; c=2.99792458e10  # cm/s
amu=1.66053907e-27; Ang=1e-10
def Qzpt_amu05Ang(wcm):
    w=2*np.pi*c*wcm            # s^-1
    q2_SI=hbar/(2*w)           # kg m^2
    # convert kg m^2 -> amu Ang^2
    q2=q2_SI/(amu*Ang**2)
    return np.sqrt(q2)         # amu^0.5 * Ang
# sanity
print("Qzpt(1000 cm-1) =",round(Qzpt_amu05Ang(1000),4),"amu^0.5*Ang  (typical ~0.1)")

eqLa=data["eq"][0]["dip"].copy()   # S1 = La at equilibrium
def La_vec(tag):
    v=data[tag][0]["dip"].copy()   # S1
    if np.dot(v,eqLa)<0: v=-v       # phase-align to eq La
    return v
def La_mu(tag):
    return data[tag][0]["mu_D"]

for m in [26,28]:
    dQ=man["modes"][str(m)]["dQ_amu05_Ang"]
    w =man["modes"][str(m)]["freq_cm"]
    vp=La_vec(f"m{m}_plus"); vm=La_vec(f"m{m}_minus")
    dmu_vec=(vp-vm)/(2*dQ)*AU2D          # D per (amu^0.5 Ang)
    dmu_mag=np.linalg.norm(dmu_vec)
    # scalar cross-check
    dmu_scalar=(La_mu(f"m{m}_plus")-La_mu(f"m{m}_minus"))/(2*dQ)
    Qz=Qzpt_amu05Ang(w)
    mu_vib=dmu_mag*Qz
    mu_vib_scalar=abs(dmu_scalar)*Qz
    print(f"\n=== Mode {m}  ({w:.1f} cm-1) ===")
    print(f"  La |mu|: eq={data['eq'][0]['mu_D']:.3f}  +={La_mu(f'm{m}_plus'):.3f}  -={La_mu(f'm{m}_minus'):.3f} D")
    print(f"  dQ(step)={dQ:.4f} amu^0.5*Ang   Qzpt={Qz:.4f} amu^0.5*Ang")
    print(f"  |dmu/dQ| (vector) = {dmu_mag:.3f} D/(amu^0.5 Ang)")
    print(f"  --> mu_vib (vector)  = {mu_vib:.3f} D")
    print(f"  --> mu_vib (scalar check) = {mu_vib_scalar:.3f} D")
    # implication for IR
    smax0=12.77; gamma=5.31; mu_e=2.1
    split=smax0*(mu_vib/mu_e)**2
    print(f"  implied dominant IR splitting @ this mu_vib = {split:.3f} cm-1   (gamma={gamma})  -> D_detect={split/gamma:.3f}")

print("\n\n########## ROBUSTEZ: mismo cálculo si La fuera S2 (Lb) ##########")
eqS2=data["eq"][1]["dip"].copy()
def S2_vec(tag):
    v=data[tag][1]["dip"].copy()
    if np.dot(v,eqS2)<0: v=-v
    return v
for m in [26,28]:
    dQ=man["modes"][str(m)]["dQ_amu05_Ang"]; w=man["modes"][str(m)]["freq_cm"]
    vp=S2_vec(f"m{m}_plus"); vm=S2_vec(f"m{m}_minus")
    dmu=np.linalg.norm((vp-vm)/(2*dQ)*2.5417464)
    Qz=Qzpt_amu05Ang(w); mv=dmu*Qz
    print(f"  Mode {m}: mu_vib(S2) = {mv:.3f} D")
