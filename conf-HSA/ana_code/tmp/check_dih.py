import math, numpy as np

def dihedral_py(p1,p2,p3,p4):
    def sub(a,b): return tuple(a[i]-b[i] for i in range(3))
    def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
    def dot(a,b): return sum(a[i]*b[i] for i in range(3))
    def norm(a): return math.sqrt(dot(a,a))
    b1=sub(p2,p1);b2=sub(p3,p2);b3=sub(p4,p3)
    n1=cross(b1,b2);n2=cross(b2,b3);m1=cross(n1,b2)
    x=dot(n1,n2)/(norm(n1)*norm(n2));y=dot(m1,n2)/(norm(m1)*norm(n2))
    return math.degrees(math.atan2(y,x))

def dihedral_np(p1,p2,p3,p4):
    b1=p2-p1;b2=p3-p2;b3=p4-p3
    n1=np.cross(b1,b2);n2=np.cross(b2,b3)
    nn1=np.linalg.norm(n1);nn2=np.linalg.norm(n2)
    n1/=nn1;n2/=nn2;b2n=b2/np.linalg.norm(b2)
    return float(np.degrees(np.arctan2(np.dot(np.cross(n1,n2),b2n),np.dot(n1,n2))))

def parse_dict(path):
    with open(path) as f: lines=f.readlines()
    natoms=int(lines[1].strip()); coords={}
    for line in lines[2:2+natoms]:
        idx=int(line[15:20])
        coords[idx]=(float(line[20:28]),float(line[28:36]),float(line[36:44]))
    return coords

def parse_np(path):
    with open(path) as f: lines=f.readlines()
    natoms=int(lines[1].strip()); xp=np.zeros((natoms,3))
    for i,line in enumerate(lines[2:2+natoms]):
        xp[i]=[float(line[20:28]),float(line[28:36]),float(line[36:44])]
    return xp

path = '/share/home/nglokwan/dparker/dp_xinyi/ana_code/invert/phe_sssD.gro'
xp_dict = parse_dict(path)
xp_np = parse_np(path)

for ai,aj,ak,al in [(19,73,74,85),(1,55,2,21)]:
    py_val = dihedral_py(xp_dict[ai],xp_dict[aj],xp_dict[ak],xp_dict[al])
    np_val = dihedral_np(xp_np[ai-1],xp_np[aj-1],xp_np[ak-1],xp_np[al-1])
    print(f"({ai},{aj},{ak},{al}): py={py_val:.3f}  np={np_val:.3f}  diff={py_val-np_val:.4f}")
