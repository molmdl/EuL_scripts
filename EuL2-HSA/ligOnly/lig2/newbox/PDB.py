import sys
import math

def create_atom(atom_id=1,
                atom_name = 'C',
                res_name = 'NON',
                chain_name = 'A',
                res_id=1,
                x=0.,y=0.,z=0.,
                seg_name = 'O'):

#    Read a line and parse it into PDB components
#    [0]: atom number
#    [1]: atom name
#    [2]: residu name
#    [3]: chain name
#    [4]: resid
#    [5-7]: x,y,z
#    [8]: segname

# how to use this routine to create a new atom in a list of atoms
#atoms.append(create_atom(atom_name= 'SI',
#                atom_id = 1,
#                res_id = 1,
#                res_name = 'PRY',
#                chain_name = 'A', x= 0.0, y=0.0, z=0.0))


 return [atom_id, atom_name, res_name, chain_name, res_id, x, y, z, seg_name]

def lookup_xp(_xp,_name,_anm):

 if not _anm in _name:
  print( 'Incorrect atom name', _anm)
  sys.exit()

 _idx = _name.index(_anm)
# print (_xp[_idx])
 return _xp[_idx]

def __nm(_s):
 return _s.split('|')[0]

def __typ(_s):
 return _s.split('|')[1]

def __typB(_s):
 return _s.split('|')[2]

def build_mol_from_pdb(fnm,res_atom,res_name='NON', res_id=1,suffix=''):
# similar to build_mol_from_txt except that
# coordinates for the frist
# three atoms will be loaded from res_atoms in pdb atom format
# first line of input will be different 
# C1|CI|C O1|OI|O C2|CA|CA
# O2|OM|O C1 C2 O1 1.34 112.1 180.0
# C2|CX|C O2 C1 C2 1.48 119.67 180.0
# here assume that C1, O1, C2 are defined in res_atom 
# which is part of a pdb that will be modified
#
# X|Y|Z X anm Y atyp Z atyp whose nbparam are used as template
# Z can be the same as Y

 _atypB = []
 _atyp=[]
 _read_first = True
 _frag = []

 _res_anm = [_i[1] for _i in res_atom]

 for rl in open(fnm, 'r'):
  srl = rl[:-1].split()
#  print(srl)
  if _read_first:
   _read_first = False
   f_3 = [[_i for _i in res_atom[_res_anm.index(__nm(srl[0]))]],\
          [_i for _i in res_atom[_res_anm.index(__nm(srl[1]))]],\
          [_i for _i in res_atom[_res_anm.index(__nm(srl[2]))]]\
         ]

#   f_3 = create_first3(__nm(srl[0]),__nm(srl[1]),float(srl[2]),\
#                      __nm(srl[3]),float(srl[4]), float(srl[5]),\
#                      res_name = res_name, res_id = res_id,\
#                      suffix = suffix)

   _atyp.append(__typ(srl[0]))
   _atyp.append(__typ(srl[1]))
   _atyp.append(__typ(srl[2]))
   _atypB.append(__typB(srl[0]))
   _atypB.append(__typB(srl[1]))
   _atypB.append(__typB(srl[2]))


   continue

  _frag.append([__nm(srl[0]),srl[1],srl[2],srl[3],\
            float(srl[4]), float(srl[5]),float(srl[6])])
  _atyp.append(__typ(srl[0]))
  _atypB.append(__typB(srl[0]))

 _atoms = create_frag(None,None,None,_frag, first_three = f_3,\
                 res_name = res_name, res_id = res_id,\
                      suffix = suffix)

 return (_atoms, _atyp, _atypB)


def build_mol_from_txt(fnm,res_name='NON', res_id=1,suffix=''):
# fnm contains info to build a mol from srcatch, an exmaple shown below
#
# C1|C O1|O 1.16 C2|CA 1.47 125.85
# O2|OM C1 C2 O1 1.34 112.1 180.0
# C2|CX O2 C1 C2 1.48 119.67 180.0
# ^
# first line: first three atoms
# other line: each for one atom
#
 _atyp=[]
 _read_first = True
 _frag = []
 for rl in open(fnm, 'r'):
  srl = rl[:-1].split()
#  print(srl)
  if _read_first:
   _read_first = False
   f_3 = create_first3(__nm(srl[0]),__nm(srl[1]),float(srl[2]),\
                      __nm(srl[3]),float(srl[4]), float(srl[5]),\
                      res_name = res_name, res_id = res_id,\
                      suffix = suffix)
   _atyp.append(__typ(srl[0]))
   _atyp.append(__typ(srl[1]))
   _atyp.append(__typ(srl[3]))
   continue
  
  _frag.append([__nm(srl[0]),srl[1],srl[2],srl[3],\
            float(srl[4]), float(srl[5]),float(srl[6])])
  _atyp.append(__typ(srl[0]))
  
 _atoms = create_frag(None,None,None,_frag, first_three = f_3,\
                 res_name = res_name, res_id = res_id,\
                      suffix = suffix)

 return (_atoms, _atyp)

def create_first3(first_nm, sec_nm,bond2, third_nm, bond3, angle3,\
       res_name='NON', res_id=1,suffix=''):
# create first three atoms with names as first_nm, sec_nm, third_nm
# atom with first_nm at (0,0,0)
# sec one at (bond2,0,0)
# third one at (sin(angle3)*bond3, cos(angle3)*bond3,0 )
  re=[]

  re.append(create_atom(atom_id = 1,\
                        atom_name = first_nm+suffix,\
                        res_name = res_name,\
                        res_id = res_id,\
                         x=0.0, y=0.0, z=0.0))

  re.append(create_atom(atom_id = 2,\
                        atom_name = sec_nm+suffix,\
                        res_name = res_name,\
                        res_id = res_id,\
                         x=bond2, y=0.0, z=0.0))

  _c = math.cos(angle3/180.0*3.141592654)
  _s = math.sin(angle3/180.0*3.141592654)
  
  re.append(create_atom(atom_id = 3,\
                        atom_name = third_nm+suffix,\
                        res_name = res_name,\
                        res_id = res_id,\
                         x= _c*bond3, y= _s*bond3, z=0.0))

  return re

def create_frag(xp1, xp2, xp3, frag, res_name='NON', res_id=1,suffix='',\
         first_three = None):
# 
# xp1,xp2,xp3: the coordinates of the first three points with which the growth of 
# the frament starts. frag contains internal geometry of the fragment in Z-matrix format.
# an example of frag is shown below:
#_frag_pry = [\
#['C1', '-1','-2','-3', 1.33, 125., 180.],\
#['C2', 'C1','-1','-2',1.33, 110., 180.],\
#['C3','C2','C1','-1',1.33, 110., 0.],\
#['C4','C3','C2','C1',1.33, 110., 0.],\
#['C5','C2','C3','C4',1.33, 120., 180.],\
#['C6','C5','C2','C1',1.33, 120., 180. ],\
#['C7','C6','C5','C2',1.33, 120.,0. ],\
#['C8','C3','C2','C1',1.33, 120.,180.],\
#['N2','C1','C2','C3',1.33, 125.,180.]
#]
# alternatively, first_three contains info for the first three atoms



 if first_three is None:
  _xp=[xp1,xp2,xp3]
  _name=['-1','-2','-3']
  re=[]
 else:
  re = [[_j for _j in _i ] for _i in first_three]
  _xp = [[_i[5], _i[6], _i[7]] for _i in re]
  _name = [_i[1] for _i in re]
 __aid = 4
 for __i in frag:
  _atom_nm = __i[0]
  _a1nm = __i[1]
  _a2nm = __i[2]
  _a3nm = __i[3]

  _x1 = lookup_xp(_xp, _name, _a1nm)
  _x2 = lookup_xp(_xp, _name, _a2nm)
  _x3 = lookup_xp(_xp, _name, _a3nm)

  _x = geom.calxyz( _x1, _x2, _x3,  __i[4],  __i[5],  __i[6])

  _xp.append(_x)
  _name.append(_atom_nm)

  re.append(create_atom(atom_id = __aid,\
                        atom_name = _atom_nm+suffix,\
                        res_name = res_name,\
                        res_id = res_id,\
                         x=_x[0], y=_x[1], z=_x[2]))
  __aid+=1

 return re


def is_idx(term):
 if term=='': return False
 for i in term:
  if not i in \
  ['0','1','2','3','4','5','6','7','8','9',' ']:
   return False
 return True

#all_system : contains all stuff, pro, sol and lip....
#sys_list: a pair [name of obj for subsystem, [list of resnames in this subsystem]]
#the residues with proper resname will be appended to each of subsystem
def seperate_system(all_system, sys_list):
 disgard_list=[]
 reslist = all_system.reslist
 for i in reslist:
  found=False
  resnm=i[2]
  
  for subsys in sys_list:
   if resnm in subsys[1]:
    found=True
    break
  if not found:
   if i[2] in disgard_list: continue
   print( i[2],all_system.atoms[i[0]][4],'is not chosen for any class of segments')
   print( 'We have following groups of residue names for different classes:')
   for kk,k in enumerate(sys_list):
    print( kk,':',k[1])
   print( 'Please choose one or this type of residue is all disgarded')
   co=sys.stdin.readline()
   if not is_idx(co[:-1]):
    print( i[2],'disgarded')
    disgard_list.append(i[2]) 
    continue
   opt=int(co)
   if opt<0 or opt>=len(sys_list): 
    print( i[2],'disgarded')
    disgard_list.append(i[2])
    continue
   subsys = sys_list[opt]
   subsys[1].append(i[2])

  for j in range(i[0],i[0]+i[1]):
   subsys[0].atoms.append(all_system.atoms[j])

 for subsys in sys_list:
  if len(subsys[0].atoms)>0:
   subsys[0].build_atomnms()
   subsys[0].res_list()
   subsys[0].build_xp()


#def show_temp_atom(x,anm):
# print_pdb_line(99999, anm, 'POPC', 'O',\
#                        3, x[0],x[1],x[2],'O66')

def print_pdb_atom(atom):
 print_pdb_line(atom[0],atom[1],atom[2],atom[3],atom[4],\
                atom[5],atom[6],atom[7],atom[8])

def pdb_atom_string(atom):
 return pdb_string(atom[0],atom[1],atom[2],atom[3],atom[4],\
                atom[5],atom[6],atom[7],atom[8])

def pdb_string(aid, anm, rnm, cnm, rid, x, y, z, snm):
 if len(anm)<3: 
  anm1=anm+' '
 else:
  anm1=anm
 if len(rnm)<4: 
  rnm1=rnm+' '
 else:
  rnm1=rnm
 if aid>99999: aid=99999
 return "ATOM  %5d%5s %4s%1s"%(aid, anm1,rnm1,cnm)+\
       '%4d    %8.3f%8.3f%8.3f'%(rid,x,y,z)+\
       '%18s%3s'%(' ',snm)

def print_pdb_line(aid, anm, rnm, cnm, rid, x, y, z, snm):
 if len(anm)<3: 
  anm1=anm+' '
 else:
  anm1=anm
 if len(rnm)<4: 
  rnm1=rnm+' '
 else:
  rnm1=rnm
 if aid>99999: aid=99999
 print( "ATOM  %5d%5s %4s%1s"%(aid, anm1,rnm1,cnm)+\
       '%4d    %8.3f%8.3f%8.3f'%(rid,x,y,z)+\
       '%18s%3s'%(' ',snm))

def n_res(atoms, st_idx, end_idx):
#return number of residues in atom idx range st_idx<=  <end_idx
 l_rid=-1
 l_nm='none'
 nres=0
 for i in range(st_idx, end_idx):
  if atoms[i][4]!=l_rid or atoms[i][2]!=l_nm:
   nres+=1
   l_rid=atoms[i][4]
   l_nm=atoms[i][2]
 return nres

def pdb_line(rl,autoAssign=False, ass_idx=0):
    """
    Read a line and parse it into PDB components
    [0]: atom number
    [1]: atom name
    [2]: residu name
    [3]: chain name
    [4]: resid
    [5-7]: x,y,z
    [8]: segname
    """
    atominfo=[]
#atom indice in PDB will not be used if autoAssign is on
#instead, all atom indices start at 0
    if (rl[0:4]!='ATOM'):
        return 'None'
    if not autoAssign:
     if not is_idx(rl[6:11]):
       atominfo.append(99999)
     else:
       atominfo.append(int(rl[6:11]))
    else:
     atominfo.append(ass_idx)

    atominfo.append(rl[12:17].strip())
    atominfo.append(rl[17:21].strip())
    atominfo.append(rl[21:22])
    atominfo.append(0 if rl[22:26]=='    ' else int(rl[22:26]))
    atominfo.append(float(rl[30:38]))
    atominfo.append(float(rl[38:46]))
    atominfo.append(float(rl[46:54]))
    if len(rl)>=75:
     kl=len(rl)
     if kl>76: kl=76 
     atominfo.append(rl[72:kl].strip())
    else:
     atominfo.append('')
    return atominfo

def load_pdb_xyz(atom):
 x=[]
 x.append(atom[5])
 x.append(atom[6])
 x.append(atom[7])
 return x

#load both atom nms and coords in certain range
def load_x_list_range(atoms, st_idx, en_idx):
 x_list=[]
 nm_list=[]
 #print st_idx, en_idx
 for i in range(st_idx, en_idx):
  #print i,atoms[i]
  x_list.append([atoms[i][5], atoms[i][6], atoms[i][7]])
  nm_list.append(atoms[i][1])
 return x_list, nm_list

#add new atom according 
#atom info in atoms
#only atom in range(st_idx, en_idx) will be checked
#return one coord
#add_hv_nmlist ['A','B','C'] three atom with name 'A', 'B' 'C'
#in range (st_idx, en_idx) of atoms will be used to
#construct the new atom
def cal_new_atom(atoms, st_idx, en_idx,\
                           add_hv_nmlist, add_hv_b, add_hv_a, add_hv_d):
 x_list, nm_list =  load_x_list_range(atoms, st_idx, en_idx)
 x_i=[]
 #print len(x_list), len(nm_list)
 for ni in add_hv_nmlist:
  x_idx= nm_list.index(ni)
  x_i.append(x_list[x_idx])
 re= geom.calxyz(x_i[0], x_i[1], x_i[2],\
            add_hv_b, add_hv_a, add_hv_d)
 return re

#add H atoms
#return list of coords
def add_new_H(atoms, st_idx, en_idx,\
                            add_h_type, add_h_nmlist):
 x_list, nm_list =  load_x_list_range(atoms, st_idx, en_idx)
 x_i=[]
 for ni in add_h_nmlist:
  x_idx= nm_list.index(ni)
  x_i.append(x_list[x_idx])
 re=[]
 if add_h_type==1:
  re.append(geom.add1H_pyra(x_i[0], x_i[1], x_i[2], x_i[3]))
 elif add_h_type==2:
  re.append(geom.calxyz(x_i[0], x_i[1], x_i[2],\
            1.0, 110.0, -120.0))
  re.append(geom.calxyz(x_i[0], x_i[1], x_i[2],\
            1.0, 110.0, 120.0))
 return re
#==========================
#load atoms including structral infos and coordinates
#from a file named as fnm
#return a list of structures, containing all infos for each atom
def load_PDB_atoms(fnm):
 atoms=[]
 f=open(fnm,"r")
 a_idx=0
 for rl in f:
  if 'END' in rl: break
  if 'ATOM' in rl or 'HETATM' in rl:
   atoms.append(pdb_line(rl, True, a_idx))
   a_idx+=1
 f.close()
 return atoms


def load_gaussian_atoms(fnm):
 atoms=[]
 f=open(fnm,'r')
 a_idx=0
 cnt=0
 for rl in f:
  if 'N-N' in rl:
   cnt+=1
   if cnt>1: 
    break
 if cnt<=1:
  print( 'Wrong Gaussian output format')
  sys.exit()
 atxt=''
 for rl in f:
  atxt+=rl[:-1].strip()
  if '@' in rl: break
 rtxt=atxt.replace('\\','|')
 st=-1
 en=-1
 cnt=0
 for i in range(len(rtxt)-1):
  if rtxt[i]=='|' and rtxt[i+1]=='|':
   cnt+=1
   if cnt==3:
    st=i+2
   if cnt==4:
    en=i

 if st<0 or en<0: 
  print( 'Wrong Guassian output format')
  sys.exit()

 ntxt=rtxt[st:en]
 srl=ntxt.split('|')
 for i in range(1,len(srl)):
  ss=srl[i].split(',')
  atoms.append([a_idx, ss[0], 'NON','X', 1, \
  float(ss[1]), float(ss[2]), float(ss[3]), 'X' ])
  a_idx+=1

 f.close()
 return atoms

class PDBdata:
 def __init__(self, fnm=None, gauss=False):
  if fnm is None:
   self.atoms=[]
   self.atomnms=[]
   self.reslist=[]
   self.xp=[]
   return
  if gauss:
   self.atoms=load_gaussian_atoms(fnm)
  else:
   self.atoms=load_PDB_atoms(fnm)
  self.build_atomnms()

  self.res_list()
  self.build_xp()

 def show_pdb(self,fnm=None,exclu=[]):
  atoms=self.atoms
  if fnm!=None:
   f=open(fnm,'w')
  for i in atoms:
   if not i[1][0] in exclu:
    if fnm!=None:
     f.write(pdb_atom_string(i)+'\n') 
    else:
     print_pdb_atom(i)
  
  if fnm!=None:
   f.write('END\n')
   f.close()
  else:
   print( 'END')

 def atoms_update(self, natoms):
  del self.reslist
  del self.xp
  del self.atomnms
  del self.atoms
# deep copy
  self.atoms=[[_j for _j in _i] for _i in natoms]
  
  self.atomnms=[]
  for i in self.atoms:
   self.atomnms.append(i[1])
  self.res_list()
  self.build_xp()


 def center(self):
  t=[0.,0.,0.]
  for x in self.xp:
   for dim,xi in enumerate(x):
    t[dim]+=  xi
  n=len(self.xp)
  for dim in range(3):
   t[dim]/=float(n)

  del self.xp
  for i in self.atoms:
   i[5]-=t[0]
   i[6]-=t[1]
   i[7]-=t[2]
  self.build_xp()

 def update_xp(self):
  xp = self.xp
  for i,j in zip(self.atoms, self.xp):
   i[5] = j[0]
   i[6] = j[1]
   i[7] = j[2]

 def update_conn(self):
  readTop.construct_conn(self, 'pdb', self.atoms)
  #print self.conn
  readTop.construct_atomSig(self)
  readTop.sort_sig(self.atomSigs)
  readTop.construct_uniq(self)


 def build_atomnms(self):
  self.atomnms=[]
  for i in self.atoms:
   self.atomnms.append(i[1])

 def build_xp(self):
  self.xp=[]
  for i in self.atoms:
   self.xp.append(load_pdb_xyz(i))

 def integrity_check(self):
  readTop.construct_conn(self,'pdb',self.atoms)
  readTop.construct_atomSig(self,None,False)
  if len(self.atomSigs[0])<len(self.atomnms)*3:
   print( 'Warning! Some part of molecule is disconnected from the other.')

#construct a residue list[[st_idx, len, name]]
 def res_list(self):
  self.reslist=[]
  atoms=self.atoms
  l_idx=-1
  for i,ai in enumerate(atoms):
   if ai[4]!=l_idx:
    if l_idx!=-1:
     self.reslist.append([st_idx, res_len, atoms[st_idx][2]])
    l_idx=ai[4]
    st_idx=i
    res_len=0
   res_len+=1
  self.reslist.append([st_idx, res_len, atoms[st_idx][2]])

#assuming xp is already defined
#this is only by reference
 def extract(self,list_nm):
  ext=[]
  for i in list_nm:
   i_idx= self.atomnms.index(i)
   ext.append(self.xp[i_idx])

  return ext

 def ext_xp_byidx(self, st, en):
  if st<0 or en>len(self.atoms):
   print( 'Wrong idx for selecting coordinates')
   sys.exit()
  xp=self.xp
  t_xp=[]
  for i in range(st, en):
   t_xp.append([xp[i][0], xp[i][1], xp[i][2]])

  return t_xp

 def extract_atompair_byres(self, list_pair):
  ext=[]
  atomnms=self.atomnms
  for i in self.reslist:
   t_nm=[]
   for j in range(i[0],i[0]+i[1]):
    t_nm.append(atomnms[j])
   for k in list_pair:
    if k[0] in t_nm and k[1] in t_nm:
     k0 = t_nm.index(k[0])+ i[0]
     k1 = t_nm.index(k[1])+ i[0]
     ext.append([[k0, self.xp[k0]],[k1, self.xp[k1]]])
  return ext

 def extract_bypart(self, list_nm, st_idx, seglen):
  ext=[]
  t_list=[self.atomnms[i] for i in range(st_idx, st_idx+seglen)]
#  print t_list
#  print list_nm

  for i in list_nm:
   if not i in t_list: continue
   i_idx=t_list.index(i) + st_idx
   ext.append(self.xp[i_idx])
  return ext

 def extractidx_bypart(self, list_nm, st_idx, seglen):
  ext=[]
  t_list=[self.atomnms[i] for i in range(st_idx, st_idx+seglen)]
#  print t_list
#  print list_nm

  for i in list_nm:
   if not i in t_list: continue
   i_idx=t_list.index(i) + st_idx
   ext.append(i_idx)
  return ext



#return number of differet chains
#and staring internal idx for each chain
 def get_chains(self,verbose=True):
  c_st=[]
  n_c=0
  l_c='$$'
  for i in self.reslist:
   st_idx=i[0]
   c_c = self.atoms[st_idx][3]
   #print c_c
   if l_c!=c_c:
    
    l_c=c_c
    n_c+=1
    if verbose:
     print( 'Chain',n_c,'starts at',\
     self.atoms[st_idx][8]+':'+self.atoms[st_idx][2]+str(self.atoms[st_idx][4]))
    c_st.append(st_idx)
  c_st.append(len(self.atoms))
  return c_st

#get both chain name and chain st_idx
#return a dict:  'chain name': [st_idx, end_idx, st_residx, end_residx]
 def get_chain_info(self):
  c_st=self.get_chains(False)
  idx_serie=[]
  chain_name=[]
  res_id=0
  for i in range(len(c_st)-1):
   e_id=res_id + n_res(self.atoms, c_st[i], c_st[i+1])
   idx_serie.append([c_st[i], c_st[i+1],res_id, e_id])
   res_id = e_id
   chain_name.append(self.atoms[c_st[i]][3])
  return dict([(i,j) for i,j in zip(chain_name, idx_serie)])
   
 def get_segments(self,verbose=True):
  c_st=[]
  n_c=0
  l_c='$$'
  for i in self.reslist:
   st_idx=i[0]
   c_c = self.atoms[st_idx][8]
   #print c_c
   if l_c!=c_c:
    
    l_c=c_c
    n_c+=1
    if verbose:
     print( 'Chain',n_c,'starts at',\
     self.atoms[st_idx][8]+':'+self.atoms[st_idx][2]+str(self.atoms[st_idx][4]))
    c_st.append(st_idx)
  c_st.append(len(self.atoms))
  return c_st

 def build_struct(self):
  if len(self.atoms)==0: return
  self.build_atomnms()
  self.res_list()
  self.build_xp()

 def tcharge(self):
  ch_res=top_db.charge_res_db()
  reslist = self.reslist
  tc=0
  for i in reslist:
   if i[2] in ch_res:
    tc+=ch_res[i[2]]
  return tc

 def volume(self):
  box = self.box()
  return box[0]*box[1]*box[2]

 def box(self):
  xp=self.xp
  max, min = geom.find_max_min(xp)
  box_c=[]
  for i,j in zip(max, min):
   box_c.append(i-j)
  return box_c

 def show_gjf(self, fnm):

  f = open(fnm,'w')

  f.write('%chk\n')
  f.write('%mem\n')
  f.write('#\n')
  f.write('\n')
  f.write('Title\n')
  f.write('\n')
  f.write('0 1\n')

  for _i in self.atoms:
   f.write('%s		%f	%f	%f\n'%\
      (_i[1][0], _i[5],_i[6],_i[7]))

  f.write('\n')
  f.write('\n')

  f.close()

 def index(self, rid, anm):
# return index(0-based) of atom with name anm in rid th (also 0-based) residue
  _rinfo = self.reslist[rid]

  _nms = self.atomnms[_rinfo[0]:_rinfo[0]+_rinfo[1]]

  if not anm in _nms: return -1

  return  _nms.index(anm) + _rinfo[0]

#end of PDBdata class

