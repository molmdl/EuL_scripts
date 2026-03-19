# Check if -d0_L (negated) == act_D for all 230 ft=2 terms
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

xp_D = parse_np('/share/home/nglokwan/dparker/dp_xinyi/ana_code/invert/phe_sssD.gro')
xp_L = parse_np('/share/home/nglokwan/dparker/dp_xinyi/ana_code/invert/phe_sssL.gro')
itp_L = open('/share/home/nglokwan/dparker/dp_xinyi/ana_code/invert/phe_sssL.itp').readlines()

STEREOCENTRE_LOCAL_SETS = [
    {74,73,85,86,87,88,92,127},{75,72,81,82,83,84,103,126},{76,71,77,78,79,80,114,125}]
SWAPPED_ATOMS = {85,127,81,126,77,125}

def is_local(ai,aj,ak,al):
    s={ai,aj,ak,al}; return any(s<=loc for loc in STEREOCENTRE_LOCAL_SETS)
def has_swapped(ai,aj,ak,al):
    return bool({ai,aj,ak,al}&SWAPPED_ATOMS)

# Test: for ALL non-local ft=2 terms, does -d0_L == act_D?
print("Testing: -d0_L == act_D for non-local ft=2 terms")
ok=0; bad=[]; swapped_ok=0; swapped_bad=[]
for line in itp_L:
    s=line.split(';')[0].strip(); parts=s.split()
    if len(parts)>=6 and parts[4]=='2':
        ai,aj,ak,al=int(parts[0]),int(parts[1]),int(parts[2]),int(parts[3])
        if is_local(ai,aj,ak,al): continue
        d0_L=float(parts[5])
        neg_d0_L=-d0_L
        # wrap
        while neg_d0_L<=-180: neg_d0_L+=360
        while neg_d0_L>180: neg_d0_L-=360
        act_D=dihedral_deg(xp_D[ai-1],xp_D[aj-1],xp_D[ak-1],xp_D[al-1])
        err=abs(neg_d0_L-act_D); err=min(err,360-err)
        if has_swapped(ai,aj,ak,al):
            if err>2: swapped_bad.append((ai,aj,ak,al,d0_L,neg_d0_L,act_D,err))
            else: swapped_ok+=1
        else:
            if err>2: bad.append((ai,aj,ak,al,d0_L,neg_d0_L,act_D,err))
            else: ok+=1

print(f"  Non-swapped: OK={ok}  BAD={len(bad)}")
print(f"  Swapped-atom: OK={swapped_ok}  BAD={len(swapped_bad)}")
if bad:
    print("  Non-swapped bad terms:")
    for r in bad[:5]: print(f"    {r}")
if swapped_bad:
    print("  Swapped bad terms:")
    for r in swapped_bad[:5]: print(f"    {r}")
