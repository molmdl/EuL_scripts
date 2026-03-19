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

xp_L = parse_np('/share/home/nglokwan/dparker/dp_xinyi/ana_code/invert/phe_sssL.gro')
itp_L = open('/share/home/nglokwan/dparker/dp_xinyi/ana_code/invert/phe_sssL.itp').readlines()

print("L-form: checking d0_L vs act_L for first 5 ft=2 terms")
print(f"{'ai':>5} {'aj':>5} {'ak':>5} {'al':>5}  {'d0_L':>9}  {'act_L':>9}  {'d0_L==act_L?':>14}  {'d0_L==-act_L?':>14}")
count=0
for line in itp_L:
    s=line.split(';')[0].strip(); parts=s.split()
    if len(parts)>=6 and parts[4]=='2':
        ai,aj,ak,al=int(parts[0]),int(parts[1]),int(parts[2]),int(parts[3])
        d0=float(parts[5])
        act=dihedral_deg(xp_L[ai-1],xp_L[aj-1],xp_L[ak-1],xp_L[al-1])
        eq = abs(d0-act)<1
        neq = abs(d0+act)<1
        print(f"{ai:5d} {aj:5d} {ak:5d} {al:5d}  {d0:9.3f}  {act:9.3f}  {'YES' if eq else 'NO':>14}  {'YES' if neq else 'NO':>14}")
        count+=1
        if count>=5: break
