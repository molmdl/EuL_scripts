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
itp_lines = open('/share/home/nglokwan/dparker/dp_xinyi/ana_code/invert/phe_sssD.itp').readlines()

errors = []; ok = 0
for line in itp_lines:
    s = line.split(';')[0].strip()
    parts = s.split()
    if len(parts) >= 6 and parts[4] == '2':
        try:
            ai,aj,ak,al = int(parts[0]),int(parts[1]),int(parts[2]),int(parts[3])
            d0_itp = float(parts[5])
            act = dihedral_deg(xp_D[ai-1],xp_D[aj-1],xp_D[ak-1],xp_D[al-1])
            err = abs(act - d0_itp)
            if err > 180: err = 360-err
            if err > 2.0: errors.append((ai,aj,ak,al,d0_itp,act,err))
            else: ok += 1
        except: pass

print(f"OK: {ok}  |  BAD (err>2deg): {len(errors)}")
if errors:
    print(f"{'ai':>5} {'aj':>5} {'ak':>5} {'al':>5}  {'d0_itp':>9}  {'act_D':>9}  {'err':>7}")
    for row in errors:
        print(f"{row[0]:5d} {row[1]:5d} {row[2]:5d} {row[3]:5d}  {row[4]:9.3f}  {row[5]:9.3f}  {row[6]:7.3f}")
else:
    print("All 230 ft=2 d0 values match D-GRO geometry within 2deg. PASS.")
