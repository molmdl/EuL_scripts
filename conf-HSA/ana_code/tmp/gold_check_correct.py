# Correct gold-standard check: d0_itp should equal act_D (not -act_D)
import numpy as np

def parse_np(path):
    with open(path) as f: lines=f.readlines()
    natoms=int(lines[1].strip()); xp=np.zeros((natoms,3))
    for i,line in enumerate(lines[2:2+natoms]):
        xp[i]=[float(line[20:28]),float(line[28:36]),float(line[36:44])]
    return xp

def dihedral_deg(p1,p2,p3,p4):
    b1=p2-p1;b2=p3-p2;b3=p4-p3
    n1=np.cross(b1,b2);n2=np.cross(b2,b3)
    nn1=np.linalg.norm(n1);nn2=np.linalg.norm(n2)
    if nn1<1e-9 or nn2<1e-9: return float('nan')
    n1/=nn1;n2/=nn2;b2n=b2/np.linalg.norm(b2)
    return float(np.degrees(np.arctan2(np.dot(np.cross(n1,n2),b2n),np.dot(n1,n2))))

# First verify L-form: d0_L should == act_L
xp_L = parse_np('/share/home/nglokwan/dparker/dp_xinyi/ana_code/invert/phe_sssL.gro')
itp_L = open('/share/home/nglokwan/dparker/dp_xinyi/ana_code/invert/phe_sssL.itp').readlines()
l_ok=0; l_bad=[]
for line in itp_L:
    s=line.split(';')[0].strip(); parts=s.split()
    if len(parts)>=6 and parts[4]=='2':
        ai,aj,ak,al=int(parts[0]),int(parts[1]),int(parts[2]),int(parts[3])
        d0=float(parts[5]); act=dihedral_deg(xp_L[ai-1],xp_L[aj-1],xp_L[ak-1],xp_L[al-1])
        err=abs(d0-act); err=min(err,360-err)
        if err>2: l_bad.append((ai,aj,ak,al,d0,act,err))
        else: l_ok+=1
print(f"L-form check: OK={l_ok}  BAD={len(l_bad)}")
if l_bad:
    for r in l_bad[:3]: print(f"  {r}")

# Now check D-form (should have d0_D = act_D, which means d0_D = -act_L = -d0_L)
xp_D = parse_np('/share/home/nglokwan/dparker/dp_xinyi/ana_code/invert/phe_sssD.gro')
itp_D = open('/share/home/nglokwan/dparker/dp_xinyi/ana_code/invert/phe_sssD.itp').readlines()
d_ok=0; d_bad=[]
for line in itp_D:
    s=line.split(';')[0].strip(); parts=s.split()
    if len(parts)>=6 and parts[4]=='2':
        ai,aj,ak,al=int(parts[0]),int(parts[1]),int(parts[2]),int(parts[3])
        d0=float(parts[5]); act=dihedral_deg(xp_D[ai-1],xp_D[aj-1],xp_D[ak-1],xp_D[al-1])
        err=abs(d0-act); err=min(err,360-err)
        if err>2: d_bad.append((ai,aj,ak,al,d0,act,err))
        else: d_ok+=1
print(f"D-form check: OK={d_ok}  BAD={len(d_bad)}")
if d_bad:
    print(f"{'ai':>5} {'aj':>5} {'ak':>5} {'al':>5}  {'d0_D':>9}  {'act_D':>9}  {'err':>7}")
    for r in d_bad:
        print(f"{r[0]:5d} {r[1]:5d} {r[2]:5d} {r[3]:5d}  {r[4]:9.3f}  {r[5]:9.3f}  {r[6]:7.3f}")
else:
    print("PASS: all 230 ft=2 d0 values match D-GRO geometry within 2deg.")
