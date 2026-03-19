import sys
import atomtyping
import readTop

#search all connection from bonding

def add_sc_list(sc, c_list):
 t=[i for i in sc]
 if t[0]>t[-1]: t.reverse()
 if not t in c_list: 
  
  c_list.append(t)
# else:
#  print 'repeat',sc

def sc_conlist(sc, mark, c12, c13, c14, direct_con):
 ne=len(sc)
 if ne==2:
  #print sc
  add_sc_list(sc, c12)
  #print c12
 if ne==3:
  add_sc_list(sc, c13)
 if ne==4:
  add_sc_list(sc, c14)
  return

 cur=sc[-1]
 for i in direct_con[cur]:
  if mark[i]==0:
   mark[i]=1
   sc.append(i)
   #print sc
   sc_conlist(sc, mark, c12, c13, c14, direct_con)
   sc.pop()
   mark[i]=0



class connect:
#restop from from Restop
 def __init__(self, fnm=None, antechamber=False,restop=None,pdb=None,psf=None):
  if fnm is None:
#reading bonding from psf obj
   if psf!=None:
    self.bond = [i[:] for i in psf.bonds]
    self.full_connect()
#reading bonding from pdb obj
   if pdb!=None:
    readTop.construct_conn(pdb,'pdb',pdb.atoms)
    
    #print self.bond
    self.full_connect(pdb.conn)
    self.bond=[[j for j in i] for i in self.c12] 
#reading bonding from residue top
   if restop!=None:
    self.bond=[] 
    for i in restop.bonds:
     #in case -X or +X for procedding or following residue
     hasnb=False
     for j in i:
      if j[0] in ['-','+']: hasnb=True
     if hasnb: continue
     t_l=[restop.atomnms.index(j) for j in i]
     self.bond.append(t_l)
    self.full_connect()
   return
#reading bonding from antechamber
  if antechamber==True:
   self.bond=[]
   f=open(fnm,'r')
   for rl in f: 
    if 'BOND' in rl:
     srl=rl[:-1].split()
     self.bond.append([int(srl[2])-1, int(srl[3])-1,\
			srl[4], srl[5]])
      
   f.close()
  self.full_connect()
  return 

 def get_direct_conn(self):
  ne=0
  for i in self.bond:
   if i[0]+1>ne: ne=i[0]+1
   if i[1]+1>ne: ne=i[1]+1
  direct_con=[[] for i in range(ne)]
  
  for i in self.bond:
   direct_con[i[0]].append(i[1])
   direct_con[i[1]].append(i[0])
  return direct_con  


 def full_connect(self, d_conn=None):
  if d_conn==None:
   direct_con = self.get_direct_conn()
  else:
   direct_con=d_conn
  ne=len(direct_con)
  
  j=0
  for i in direct_con: j+=len(i)
  #print j

  c12=[]
  c13=[]
  c14=[]
  #print ne
  for i in range(ne):
   mark = [0 for j in range(ne)]
   mark[i]=1
   sc=[i]
   sc_conlist(sc, mark, c12, c13, c14, direct_con)

  self.c12=c12
  self.c13=c13
  self.c14=c14

#genertae impr list according to hybridization states
 def gen_impr_list(self, atomnms,amber=False):
  t=[]
  conn = self.get_direct_conn()
  hyd_state = atomtyping.get_hybrid_state(atomnms, conn)
  for i,j in enumerate(hyd_state):
   if j in ['Nsp2','Csp2'] and len(conn[i])==3:
    if amber:
     c_t=[]
     for k in conn[i]: c_t.append(k)
     c_t.insert(2,i)
    else:
     c_t=[i]
     for k in conn[i]:
      c_t.append(k)
    t.append(c_t)

  self.impr = t

def crowl_molecule(conn, atomnms): #return the list of atom order following the order as crowled
 q=[]
 head=0
 q.append(0)
 natom=len(atomnms)
 mark=[0 for i in range(natom)]
 mark[0]=1
 crl_mol(q, 0, mark, conn, atomnms)

 return q


def crl_mol(q, idx, mark, conn, atomnms):
 for i in conn[idx]:
  if mark[i]==0 and atomnms[i][0]=='H':
   q.append(i)
   mark[i]=1
 for i in conn[idx]:
  if mark[i]==0:
   q.append(i)
   mark[i]=1
   crl_mol(q, i, mark, conn, atomnms)

