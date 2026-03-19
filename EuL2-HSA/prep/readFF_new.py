from readTop import *
import sys
import top_db_new
from out_put import *
import geom
import helper 
from atomtyping import *
import readTop
import numpy as np
import geom_opt

def is_val_diff(vA,vB ,rel_tor):

 if vA==0.0 and vB==0.0: return False
 

 _b = max(abs(vA), abs(vB))

 _d = abs(vA-vB)

 if _d/_b> rel_tor:
  return True
 else:
  return False




def est_nb_matr(ff):
#populated nonbonded parameter matrix for ff
#return _is_nb: matrix indicating if there are nb parameters between types i and j
#       _c6_nb: C6 param for i-j
#       _c12_nb: C12 param for i-j

 _atom = [_i[0] for _i in ff.gmx_atoms]
 _n = len(_atom)

 _is_nb = np.zeros((_n, _n), dtype=np.int32)
 _c6_nb = np.zeros((_n, _n))
 _c12_nb = np.zeros((_n, _n))

 for _i, _j in zip(ff.gmx_nbs, ff.gmx_nb_params):

  _a = _i[0]
  _b = _i[1]

  if not _a in _atom:
   print(_a, 'not defined ***')
   continue
  if not _b in _atom:
   print(_b, 'not defined ***')
   continue

  _idx_a = _atom.index(_a)
  _idx_b = _atom.index(_b)

  _is_nb[_idx_a, _idx_b] = 1
  _is_nb[_idx_b, _idx_a] = 1
 
  _c6_nb[_idx_a, _idx_b] = _j[0]
  _c6_nb[_idx_b, _idx_a] = _j[0]

  _c12_nb[_idx_a, _idx_b] = _j[1]
  _c12_nb[_idx_b, _idx_a] = _j[1]

 return (_is_nb, _c6_nb, _c12_nb)

def compare_nbparam(ffA, ffB, compare_list=None, rel_tor=1e-4):
# compare the nonbonded parameters between two force fields ffA and ffB
# in two modes: 1, if compare_list is specified with a list of atom types
# then the parameters associated with these types will be compared between the two forcefields
# 2, all the atom types in ffA will be used for comparison instead

 if compare_list is None:
   print('Comparison between two force fields for types of the first force field')
   c_list = [_i[0] for _i in ffA.gmx_atoms]
 else:
   c_list = compare_list

 _nc = len(c_list)

 _atomA_list = [_i[0] for _i in ffA.gmx_atoms]
 _atomB_list = [_i[0] for _i in ffB.gmx_atoms]

 print('checking force field A...')
 _is_nbA, _c6_nbA, _c12_nbA = est_nb_matr(ffA)
 print('checking force field B...')
 _is_nbB, _c6_nbB, _c12_nbB = est_nb_matr(ffB)
 
 for _i in range(_nc):
  if not c_list[_i] in _atomA_list:
    print(c_list[_i],'not found in force field A')
    continue
  if not c_list[_i] in _atomB_list:
    print(c_list[_i],'not found in force field B')
    continue

  _idx_A_i = _atomA_list.index(c_list[_i])
  _idx_B_i = _atomB_list.index(c_list[_i])

  for _j in range(_i, _nc, 1):
   if not c_list[_j] in _atomA_list:
    continue
   if not c_list[_j] in _atomB_list:
    continue

    
   _idx_A_j = _atomA_list.index(c_list[_j])
   _idx_B_j = _atomB_list.index(c_list[_j])

   if _is_nbA[_idx_A_i, _idx_A_j]==0:

    if _is_nbB[_idx_B_i, _idx_B_j]==0:
     continue
    else:
     print('nb param for',c_list[_i],c_list[_j],'defined in force field B but not in A')
     continue
   elif _is_nbB[_idx_B_i, _idx_B_j]==0:
    print('nb param for',c_list[_i],c_list[_j],'defined in force field A but not in B')
    continue

   _c6A = _c6_nbA[_idx_A_i, _idx_A_j]
   _c12A = _c12_nbA[_idx_A_i, _idx_A_j]
   _c6B = _c6_nbB[_idx_B_i, _idx_B_j]
   _c12B = _c12_nbB[_idx_B_i, _idx_B_j]

   if is_val_diff(_c6A, _c6B,rel_tor):
    print('C6 param for',c_list[_i],c_list[_j],'in A',_c6A,': in B',_c6B)

   if is_val_diff(_c12A, _c12B,rel_tor):
    print('C12 param for',c_list[_i],c_list[_j],'in A',_c12A,': in B',_c12B)




#convert parameters according to rules
def converParam(pi, p_field, convert):
 re=[]
 for pidx,i,j in zip(range(len(p_field)), p_field, convert):
  
  if i==-1:
   re.append(j)
  elif i==-2:
   re.append(pi[pidx])
  else:
   re.append(pi[i]*j)
 return re

#wild card will not considered
def is_dihe_in_FF(ityp, multi, dihes, dihe_params):
 a=[i for i in ityp]
 ra=[i for i in ityp]
 ra.reverse()
 for i,j in zip(dihes, dihe_params):
  if a==i or ra==i:
   if multi==j[1]:
    return True
 return False


def idx_dihe_in_FF(ityp, multi, dihes, dihe_params):
 a=[i for i in ityp]
 ra=[i for i in ityp]
 ra.reverse()
 idx=0
 for i,j in zip(dihes, dihe_params):
  if a==i or ra==i:
   if multi==j[1]:
    return idx
  idx+=1
 return -1



#generate of idx entry of impr for l
#l is atom idx multiplet not typ
#need to be translated with t_typ
#note: l also will be changed according to best match with FF
def match_impr_amber_list(l, imprs, t_typ):
 t=[]
 of_idx=[[0,1,3],[0,3,1],\
         [1,0,3],[1,3,0],\
         [3,0,1],[3,1,0]]
 found=0
 for fi in of_idx:
  c_t=[t_typ[l[fi[1]]],\
       t_typ[l[fi[0]]],\
       t_typ[l[2]],\
       t_typ[l[fi[2]]]]
  for j,k in enumerate(imprs):
   if match_entry(c_t, k):
    t.append(j)
    n_l = [l[fi[1]],\
           l[fi[0]],\
           l[2],\
           l[fi[2]]]
    l=n_l
    found=1
    break
  if found==1:
   break
 return t
  
def match_entry(l,f):
 if len(l)!=len(f): return False
 for i,j in zip(l,f):
  if j=='X': continue
  if i!=j: return False
 return True

def match_FF_entry(l,f): # return True if l matched f, a FF entry
 c_l=[i for i in l]
 if match_entry(c_l, f):
  return True
 c_l.reverse()
 if match_entry(c_l, f):
  return True
 return False
#match with *
def match_FF_list_wildcard(t_l, flist):
 t=[]
 k_l=[i for i in t_l]
 k_l.reverse()
 for ii,i in enumerate(flist):
  a_l=[]
  b_l=[]
  c_l=[]
  d_l=[]
  for si, sk, sj,sl in zip(t_l,k_l, i,i):
   if '*' in si:
    idx=si.index('*')
    hi=si[:idx]
    si=hi
    if sj!='X':
     hj=sj[:idx]
     sj=hj

   if '*' in sk:
    idx=sk.index('*')
    hk=sk[:idx]
    sk=hk
    if sl!='X':
     hl=sl[:idx]
     sl=hl


   a_l.append(si)
   b_l.append(sj)
   c_l.append(sk)
   d_l.append(sl)
   
  if match_entry(a_l, b_l) or match_entry(c_l,d_l):
   t.append(ii)
 return t

def match_FF_list(l, flist): #return list of FF entry idx for matched multiplet
 t=[]
 
 for i,f in enumerate(flist):
  if match_FF_entry(l, f): t.append(i)
 return t
 
def sort_reverse(l):
 if l[0]>l[-1]: l.reverse()

def sort_reverse_list(ls):
 for i in ls:
  sort_reverse(i)


def load_mass(fnm):

 _atyp = []
 _mass = []

 for _rl in open(fnm, 'r'):

  srl =_rl[:-1].split()
  if len(srl)==0:continue
  if srl[0]=='MASS':
   _atyp.append(srl[2])
   _mass.append(float(srl[3]))

 return dict([(_i,_j) for _i,_j in zip(_atyp, _mass)])

class FFdata:
#read parameters either from file
#or topologies
#SCNB when used to scale 14 pair
#disgard ids of col disgarded when reading atom section in GMX format
#isdirect load nbfix parameters directly
 def __init__(self, fnm=None, charmm=True, gmxtops=None, gmxFF=False, SCNB=1.0, FF=None,disgard=[],isdirect=False):

  self.bonds=[]
  self.bond_params=[]
  self.angles=[]
  self.angle_params=[]
  self.dihes=[]
  self.dihe_params=[]
  self.imprs=[]
  self.impr_params=[]
  self.nbs=[]
  self.nb_params=[]
  self.nbfixs=[]
  self.nbfix_params=[]
  self.header=[]
  self.mass=[]
  if fnm is None: return

  if not gmxtops is None:
   self.add_param_from_gmxtop(gmxtops)
  elif gmxFF:
   self.get_param_from_gmx_file(fnm, SCNB, FF, disgard, isdirect, charmm=charmm)
  elif charmm:
   self.get_param_from_charmm_file(fnm)
   self.mass = load_mass(fnm) # return a dictionary of mass as atyp
  else:
   print( 'Err! Cannot load parameters')
   sys.exit()

  self.redundancy_check()

  

#====end of init
 def generate_nbfix(self):
  for i in self.nbfix_params:
   emin = -i[0]
   rmin = i[1]
   emin14 = -i[2]
   rmin14 = i[3]
   C6, C12 = EminRmin2Coff(emin, rmin)
   C6_14, C12_14 = EminRmin2Coff(emin14, rmin14)
   i[0]= C12
   i[1]= C6
   i[2]= C12_14
   i[3]= C6_14 
  return

 def redundancy_check(self):
  return

#loading parameters, we can use other force field FF, must be in gmx format
#as a ref to generate pair parameters
#note this function append new parameters
 def get_param_from_gmx_file(self, fnm, SCNB, FF=None,disgard=[],isdirect=False,\
#         charmm force field format in GMX 
           charmm=False):
  print('col',disgard,'will be ignored from',fnm)
  self.gmx_atoms=[]
  self.gmx_nbs=[]
  self.gmx_pairs=[]
  self.gmx_dihes=[]
  self.gmx_angles=[]
  self.gmx_bonds=[]
  self.gmx_atom_params=[]
  self.gmx_nb_params=[]
  self.gmx_pair_params=[]
  self.gmx_dihe_params=[]
  self.gmx_angle_params=[]
  self.gmx_bond_params=[]
  self.gmx_urdihes=[]
  self.gmx_urdihe_params=[] 
  self.gmx_imprs=[]
  self.gmx_impr_params=[] 
 
  load_rule = top_db_new.load_param_gmx_rule(self, charmm=charmm)
  all_modes = load_rule.keys()
  sel_modes = [__i for __i in all_modes]
  mute_modes = ['constrainttypes']
  mode=None
  f=open(fnm,'r')
  for _rl in f:
   rl = _rl.replace(';',' ;')
   srl=rl[:-1].split()
   if len(srl)==0: continue
   if srl[0][0]=='#': continue
   if srl[0][0]==';':
    #self.header.append(rl[:-1])
    continue
   
   mode_t = readTop.word_in_line(rl, sel_modes+mute_modes)
   if not mode_t is None:
    mode = mode_t
    continue
   if mode in mute_modes: continue

   if not mode is None:
    if mode=='atomtypes':
     srl1=[]
     for di,dj in enumerate(srl):
      if not di in disgard:
       srl1.append(dj)
     srl=srl1 
#    print(fnm,mode,srl)
    load_gmx_line(srl, mode, load_rule)
  f.close()
  #here convert gmx param to charmm

  #rule_key=load_rule.keys()
  #for ik in rule_key:
  # print ik
  # for i,j in zip(load_rule[ik][0], load_rule[ik][2]): print i,j
  sec_list=['atomtypes','nonbond-params','pairtypes',\
         'dihedraltypes','URdihe','impropertypes','bondtypes','angletypes']
  chr_list=['NONBONDED','NBFIX','NBFIX','DIHEDRALS',\
              'DIHEDRALS','IMPROPER','BONDS','ANGLES']
  load_rule_charmm = top_db_new.load_param_charmm_rule(self)
  if not FF is None:
   typnmFF=[i[0] for i in FF.nbs]
  else:
   typnmFF=[]
  #print typnmFF
  typnm=[i[0] for i in self.nbs]
#  nb_list = []
  for isl,csl in zip(sec_list, chr_list):
  # if isl in load_rule:
   l_r=load_rule[isl] 
   ftyp_gmx = l_r[0]
   fparam_gmx = l_r[2]
   l_r_chr = load_rule_charmm[csl]
   ftyp = l_r_chr[0]
   fparam = l_r_chr[2]
   for tp,param in zip(ftyp_gmx, fparam_gmx):

    if isl=='atomtypes' and (not isdirect):
     ftyp.append(tp[:])
     if charmm:
      sigma =param[2]
      epsilon = param[3]
      C6 = 4.0*epsilon*sigma**6.0
      C12 = 4.0*epsilon*sigma**12.0
      emin, rmin = geom_opt.Coff2EminRmin(C6,C12)
     else:
      emin, rmin = geom_opt.Coff2EminRmin(param[2],param[3])
     emin /= 4.184
     rmin *= 10.
     fparam.append([0.0,-emin, rmin/2. ])
     typnm.append(tp[0])
     #print typnm

    if isl=='atomtypes' and isdirect:
     ftyp.append(tp[:])
     #emin, rmin = geom.Coff2EminRmin(param[2],param[3])
     #emin /= 4.184
     #rmin *= 10.
     C12 = param[3]* 1.e+12/ 4.184
     C6  = param[2]* 1.e+6/ 4.184
     fparam.append([0.0, C12, C6 ])
     typnm.append(tp[0])


    if isl=='nonbond-params':
# or isl=='pairtypes':
     tp.sort() #since their would be only two elements in tp
     if tp in ftyp:
      print('Warning: duplicate',isl,csl,' parameters for',tp,'will be ignored.')
      continue
     else:
      ftyp.append(tp[:])

     if charmm:
      sigma = param[0]
      epsilon = param[1]
      C6 = 4.0*epsilon*sigma**6.0
      C12 = 4.0*epsilon*sigma**12.0
      emin, rmin = geom_opt.Coff2EminRmin(C6,C12)
     else:
      emin, rmin = geom_opt.Coff2EminRmin(param[0],param[1])

     emin /= 4.184
     rmin *= 10.

     fparam.append([-emin, rmin])


    if isl=='pairtypes':
# or isl=='pairtypes':
     tp.sort() #since their would be only two elements in tp
     if tp in ftyp:
      _idx_tp = ftyp.index(tp)
     else:
      ftyp.append(tp[:])
      #print tp[0], typnm
      idxFFa=-1
      idxFFb=-1
      if tp[0] in typnm: 
       idxa=typnm.index(tp[0])
      elif tp[0] in typnmFF:
       idxFFa=typnmFF.index(tp[0])
      else:
       print( 'No LJ parameters for',tp[0])
       sys.exit()
      if tp[1] in typnm: 
       idxb=typnm.index(tp[1])
      elif tp[1] in typnmFF:
       idxFFb=typnmFF.index(tp[1])
      else:
       print( 'No LJ parameters for',tp[1])
       sys.exit()

      if idxFFa>=0:
       pa = FF.nb_params[idxFFa]
      else:
       pa = self.nb_params[idxa]
    
      if idxFFb>=0:
       pb = FF.nb_params[idxFFb]
      else:
       pb = self.nb_params[idxb]


     # if not isdirect:
      emin = (pa[1]*pb[1])**0.5
      rmin = pa[2]+pb[2]
      fparam.append([-emin, rmin])
      _idx_tp = ftyp.index(tp)

     _p = fparam[_idx_tp]
     if len(_p)==4:
#      print('Warning: duplicate',isl,csl,' parameters for',tp,'will be ignored.')
      continue
     if len(_p)<2:
      print('Error: something wrong with',isl,csl,'paraneters for',tp,':',_p)
      sys.exit()
     
     if charmm:
      sigma = param[0]
      epsilon = param[1]
      C6 = 4.0*epsilon*sigma**6.0
      C12 = 4.0*epsilon*sigma**12.0
      emin14, rmin14 = geom_opt.Coff2EminRmin(C6,C12)
     else:
      emin14, rmin14 = geom_opt.Coff2EminRmin(param[0],param[1])

     emin14 /= 4.184
     rmin14 *= 10.

     _p.append(-emin14)
     _p.append(rmin14)


    #  else:
       #emin, rmin = geom.Coff2EminRmin(param[2],param[3])  
    #   emina, rmina = geom.Coff2EminRmin(pa[2],pa[1])
    #   eminb, rminb = geom.Coff2EminRmin(pb[2],pb[1])
    #   emin = (emina*eminb)**0.5
    #   rmin = (rmina+rminb)/2.
    #   fparam.append([4*emin*rmin**12,4*emin*rmin**6,\
    #                  4*emin14*rmin14**12,4*emin14*rmin14**6])

#     if isl=='nonbond-params' and (not isdirect):
#      emin, rmin = geom.Coff2EminRmin(param[0], param[1])
#      emin /= 4.184
#      rmin *= 10.
#      idx = ftyp.index(tp)
#      fparam[idx][0] = -emin
#      fparam[idx][1] = rmin
#     elif isl=='pairtypes' and (not isdirect):
#      emin, rmin = geom.Coff2EminRmin(param[0], param[1])
#      emin /= 4.184
#      rmin *= 10.
#      idx = ftyp.index(tp)
#      fparam[idx][2] = -emin
#      fparam[idx][3] = rmin
#     elif isl=='nonbond-params' and isdirect:
#      idx = ftyp.index(tp)
#      fparam[idx][0] = param[1]*1.e+12/4.184
#      fparam[idx][1] = param[0]*1.e+6/ 4.184
#     elif isl=='pairtypes' and  isdirect:
#      idx = ftyp.index(tp)
#      fparam[idx][2] = param[1]*1.e+12/4.184
#      fparam[idx][3] = param[0]*1.e+6/ 4.184


    if isl=='bondtypes':
     ftyp.append(tp[:])
     fparam.append([param[1]/2./100./4.184, param[0]*10.])
 
    if isl=='angletypes':
     ftyp.append(tp[:])
     if charmm:
      fparam.append([param[1]/2./4.184, param[0], param[3]/2./4.184, param[2]*10.0])
     else:
      fparam.append([param[1]/2./4.184, param[0]])

    if isl=='dihedraltypes':
     if charmm:
      ftyp.append(tp[:])
     else:
      ftyp.append(['X']+tp[:]+['X'])
     fparam.append([param[1]/4.184, param[2], param[0]])
# need to add impr

    if isl=="impropertypes":
     ftyp.append(tp[:])
     fparam.append([param[1]/2./4.184,0, param[0]])

    if isl=='URdihe':
     KC = geom.ur2dihe(param)
     for pi, ki in zip(range(1,6), KC):
      ftyp.append(['X']+tp[:]+['X'])
      fparam.append([ki/4.184, pi, 0.0])



 def append_gmx(self, ff):
  self.gmx_atoms+= ff.gmx_atoms
  self.gmx_nbs+= ff.gmx_nbs
  self.gmx_pairs+= ff.gmx_pairs
  self.gmx_dihes+= ff.gmx_dihes
  self.gmx_angles+= ff.gmx_angles
  self.gmx_bonds+= ff.gmx_bonds
  self.gmx_atom_params+= ff.gmx_atom_params
  self.gmx_nb_params+= ff.gmx_nb_params
  self.gmx_pair_params+= ff.gmx_pair_params
  self.gmx_dihe_params+= ff.gmx_dihe_params
  self.gmx_angle_params+= ff.gmx_angle_params
  self.gmx_bond_params+= ff.gmx_bond_params
  self.gmx_urdihes+= ff.gmx_urdihes
  self.gmx_urdihe_params+= ff.gmx_urdihe_params

 

 def get_param_from_charmm_file(self, fnm):
  load_rule = top_db_new.load_param_charmm_rule(self)
  all_modes = load_rule.keys()
  mode=None
  f=open(fnm,'r')
  for rl in f:
   srl=rl[:-1].split()
   if len(srl)==0: continue
   if srl[0][0]=='!': continue
   if 'cutnb' in srl: continue
   if 'HBOND' in srl: continue
   if 'END' in srl: continue
   if 'CTOFNB' in srl: continue
   if srl[0][0]=='*':
    self.header.append(rl[:-1])
    continue
   mode_t = readTop.word_in_line(rl, all_modes)
   if not mode_t is None:
    mode = mode_t
    continue
   if not mode is None:
    load_charmm_line(srl, mode, load_rule)

  f.close() 
#appending FF
 def append(self, FF):
  load_rule_self = top_db_new.load_param_charmm_rule(self) 
  load_rule_FF = top_db_new.load_param_charmm_rule(FF)
  tkeys = load_rule_FF.keys()
  for ik in tkeys:
   if ik=='CMAP': continue
   l_r_self = load_rule_self[ik]
   l_r_FF = load_rule_FF[ik]

   ftyp_self= l_r_self[0]
   ftyp_FF=l_r_FF[0]
   fparam_self=l_r_self[2]
   fparam_FF = l_r_FF[2]
   for itp, ipar in zip(ftyp_FF, fparam_FF):
    a=itp[:]
    ra=itp[:]
    ra.reverse()
    ok=True
    if a in ftyp_self or ra in ftyp_self:
     if ik=='DIHEDRALS':
      if is_dihe_in_FF(a, ipar[1], ftyp_self, fparam_self):
       ok=False
     else: 
      ok=False
    if not ok:
     print( ik,a,'duplicated! Disgarded')
    else:
     ftyp_self.append(itp[:])
     fparam_self.append(ipar[:])
 
#====
 def add_param_from_gmxtop(self,gmxtops):
  for itop in gmxtops:
   self.add_param_from_gmx_top(itop)


 def add_param_from_gmx_top(self,itop):
  adapt_top_2_param = top_db_new.gmx_top_param_connection(itop)
  plug_param = top_db_new.load_param_charmm_rule(self)

  
  opts_top = adapt_top_2_param.keys()
  attr = itop.atom_attr
  n_2_t_lkup =dict([(i, attr[i][0]) for i in attr.keys() ])
  for iopt in opts_top:
   if not iopt in plug_param: 
    print( 'Warning! In',iopt,'not found in the forcefield object')
    return
   self.add_param_gmx_one_sec(adapt_top_2_param[iopt],\
   plug_param[iopt],n_2_t_lkup, iopt)


 def add_param_gmx_one_sec(self, i_adapt, i_plug, lkup,iopt='None'):
  top_nm = i_adapt[0]
  top_p =i_adapt[1]
  ff_nm=i_plug[0]
  ff_p = i_plug[2]
  map_rule = i_adapt[2]
  convert = i_adapt[3]
  typ_rule = i_plug[4]
  for i_nm,i_p in zip(top_nm, top_p):
   #print i_nm, i_p
   ok=True
   if i_p==[]: continue   
   il_nm = [lkup[ii] for ii in i_nm]
   ir_nm = il_nm[:]
   ir_nm.reverse()
   rep_idx=-1
   if il_nm in ff_nm or ir_nm in ff_nm:
    if iopt=='DIHEDRALS':
     if is_dihe_in_FF(il_nm, int(i_p[2]), ff_nm, ff_p):
      rep_idx = idx_dihe_in_FF(il_nm, int(i_p[2]), ff_nm, ff_p)
       
    else:
     if il_nm in ff_nm: 
      rep_idx=ff_nm.index(il_nm)
     else:
      rep_idx=ff_nm.index(ir_nm)
     
    #if not ok:
     #for kk,ll in  zip(ff_nm, ff_p): print kk,ll
     #print iopt, int(i_p[2])
     #print 'Warning! multiplet ',i_nm, 'has duplicated parameters with ', il_nm,'in FF'
     #continue 
   # ok=False
    #continue

   t_l=[]
   for i,j,t in zip(map_rule, convert, typ_rule):
    if i<0:
     z=j
    else:
     if t==0:
      r=int(i_p[i])
     else:
      r=float(i_p[i])
     z= r*j
    t_l.append(z)
   if rep_idx>=0:
    if not helper.is_float_list_same(t_l, ff_p[rep_idx]):
     ok=False
   if ok and rep_idx<0:
    ff_nm.append(il_nm)
    ff_p.append(t_l)
   elif not ok:
    print( 'warning!',iopt,i_nm,'has specified parameters',t_l,'for',il_nm)
    print( 'In forcefield, they are different:',ff_p[rep_idx],'Disgard the specified parameters')
    #print 'overlap', il_nm,t_l
    pass

#===

 def show(self,title, fnm=None, isdirect=False):
# warning: output for nb and angle is not correct
  load_rule = top_db_new.load_param_charmm_rule(self)
  to_output=['BONDS','ANGLES','DIHEDRALS','IMPROPER',\
             'NONBONDED','NBFIX']

  if not fnm is None:
   f=open(fnm, 'w')
  else:
   f=None
 
  print_f(f,'*'+title)
  for i in self.header:
   print_f(f, i)
  print_f(f)
 
  for iout in to_output:
   if iout in load_rule:
    print_f(f, iout)
    l_r = load_rule[iout]
    param_col = l_r[3]
    param_typ = l_r[4]
    ff_nm = l_r[0]
    ff_param = l_r[2]
    for inm, iparam in zip(ff_nm, ff_param):
     p_str=''
     for i in inm: p_str+= '%6s '%(i)
     for i,j in zip(iparam, param_typ):
      if j==0:
       p_str+=str( i)+' '
      elif j==1:
       p_str+= '%10.6g '%(i)
     if iout=='NBFIX' and isdirect:
      p_str+='   direct'
     print_f(f, p_str)
   print_f(f)
  if not fnm is None:
   f.close()
#===

 def show_in_gmx(self,fnm, dihtyp = 9, angtyp= 5,\
                LJrule =2):
# for charm forcefield, in GMX, its dih type is 9
# its ang type is 5
  
  _f=open(fnm,'w')

  _f.write('[ atomtypes ]\n')
  _all_atyps= []
  _param_atyps=[]
  for _i,_j in zip(self.nbs, self.nb_params):
   _atyp = _i[0]

   C6, C12 = geom_opt.EminRmin2Coff(-_j[1]*4.184, _j[2]*0.1*2)
   sigma, epsilon = geom_opt.Coff2sigmaeps(C6, C12)

   if not _atyp in _all_atyps:
    _all_atyps.append(_atyp)
    if len(_j)==3:
     _param_atyps.append([_j[1],_j[2]])
    else:
     _param_atyps.append([_j[4],_j[5]])
   else:
    print('Warning: duplicate def of atyp',_atyp)
    print('Warning: its parameters e, rmin', _j[1],_j[2],'will be ignored.')
    continue

   if LJrule==2:
    _f.write('%s  0  %f  0.0  A  %f  %f\n'%\
     (_atyp, self.mass[_atyp], sigma, epsilon))
   elif LJrule==1:
    _f.write('%s  0  %f  0.0  A  %g  %g\n'%\
     (_atyp, self.mass[_atyp], C6, C12))


  _f.write('[ bondtypes ]\n')
  for _i,_j in zip(self.bonds, self.bond_params):
#   if 'X' in _i:
#    _f.write('; ')
   _f.write('%s  %s  1  %f  %f\n'%\
     (_i[0],_i[1], _j[1]*0.1, _j[0]*2.0*100.0*4.184))

  _f.write('[ angletypes ]\n')
  for _i,_j in zip(self.angles, self.angle_params):
#   if 'X' in _i:
#    _f.write('; ')
   if len(_j)==2:
    _f.write('%s  %s  %s  %d  %f  %f  0.0  0.0\n'%\
     (_i[0],_i[1],_i[2], angtyp , _j[1], _j[0]*2.0*4.184))
   else:
    _f.write('%s  %s  %s  %d  %f  %f  %f  %f\n'%\
     (_i[0],_i[1],_i[2], angtyp , _j[1], _j[0]*2.0*4.184,\
      _j[3]*0.1, _j[2]*2.0*4.184 ))



  _f.write('[ dihedraltypes ]\n')
  for _i,_j in zip(self.dihes, self.dihe_params):
#   if 'X' in _i:
#    _f.write('; ')
   _f.write('%s  %s  %s  %s  %d  %f  %f  %d\n'%\
     (_i[0],_i[1],_i[2],_i[3],dihtyp,\
            _j[2], _j[0]*4.184, _j[1]))

  _f.write('[ dihedraltypes ]\n')
  for _i,_j in zip(self.imprs, self.impr_params):
#   if 'X' in _i:
#    _f.write('; ')
   _f.write('%s  %s  %s  %s  2  %f  %f\n'%\
     (_i[0],_i[1],_i[2],_i[3],\
            _j[2], _j[0]*2.0*4.184))

# nonbonfix
  _f.write('[ nonbond-params ]\n')
  for _i,_j in zip(self.nbfixs, self.nbfix_params):
   C6, C12 = geom_opt.EminRmin2Coff(-_j[0]*4.184, _j[1]*0.1)
   sigma, epsilon = geom_opt.Coff2sigmaeps(C6, C12)

   if LJrule==2:
    _f.write('%s  %s  1  %f  %f\n'%(_i[0],_i[1],sigma, epsilon))
   elif LJrule==1:
    _f.write('%s  %s  1  %g  %g\n'%(_i[0],_i[1], C6, C12))

# pairs
  _f.write('[ pairtypes ]\n')
#  _f.write(';from combination of sedonadry nonbond LJ\n')
  n_atyp=len(_all_atyps)
  _e = np.zeros((n_atyp, n_atyp))
  _rmin = np.zeros((n_atyp, n_atyp))
  for _i,_j in zip(self.nbs, self.nb_params):
   _atyp = _i[0]
   _idx = _all_atyps.index(_atyp)
   if len(_j)==3: continue   
   _e_i = _j[4]
   _rmin_i = _j[5]
   for _k in range(n_atyp):
    _e_k = _param_atyps[_k][0]
    _rmin_k = _param_atyps[_k][1]
    _e_comb = - (_e_i*_e_k)**0.5
    _rmin_comb = _rmin_i + _rmin_k
    _e[_idx,_k] = _e_comb
    _rmin[_idx,_k] = _rmin_comb
    _e[_k,_idx] = _e_comb
    _rmin[_k,_idx] = _rmin_comb

# check special pair from NBFIX
# will override those from combination of LJ parameters
  for _i,_j in zip(self.nbfixs, self.nbfix_params):
   if len(_j)==2: continue
   _idx_i = _all_atyps.index(_i[0])
   _idx_j = _all_atyps.index(_i[1])
   _e_ij = _j[2]
   _rmin_ij = _j[3]
   _e[_idx_i, _idx_j] = _e_ij
   _rmin[_idx_i,_idx_j] = _rmin_ij
   _e[_idx_j, _idx_i] = _e_ij
   _rmin[_idx_j,_idx_i] = _rmin_ij

  for _i in range(n_atyp):
   for _j in range(_i, n_atyp,1):
    if _e[_i,_j]!=0.0 or _rmin[_i,_j]!=0.0:
      C6, C12 = geom_opt.EminRmin2Coff(-_e[_i,_j]*4.184, _rmin[_i,_j]*0.1)
      sigma, epsilon = geom_opt.Coff2sigmaeps(C6, C12)
      if LJrule==2:
       _f.write('%s  %s  1  %f  %f\n'%(_all_atyps[_i],_all_atyps[_j],sigma, epsilon))
      elif LJrule==1:
       _f.write('%s  %s  1  %g  %g\n'%(_all_atyps[_i],_all_atyps[_j], C6, C12))


  _f.close()

#==end of class FFdata


def load_gmx_line(srl, mode, load_rule):
 #print mode, srl
 l_r= load_rule[mode]
 l_type_pos= l_r[5]
 l_sub=l_r[6]
 if l_type_pos >= helper.len_comment(srl):
  print( 'Err! wrong line for',mode,'at',srl)
  sys.exit()
 ta = srl[l_type_pos]
 if ta in l_sub:
  tb=l_sub[ta]
  load_charmm_line(srl, tb, load_rule)
  return
 load_charmm_line(srl, mode, load_rule)

#load charmm parameters from a text line
def load_charmm_line(srl, mode, load_rule):
 l_r = load_rule[mode]
# for now, CMAP is not supported!
 if mode=='CMAP': return
 other_op = l_r[5]
 param_col = l_r[3]
 param_typ = l_r[4]
 ff_nm = l_r[0]
 ff_nm_len = l_r[1]
 ff_param = l_r[2]
 #print srl
 if helper.len_comment(srl)<=max(param_col):
  _other_op_fine = False
  if other_op:
   param_col = l_r[6]
   param_typ = l_r[7]
   _other_op_fine= True
   if helper.len_comment(srl)<=max(param_col):
    _other_op_fine = False
  if not _other_op_fine: 
   print( 'Err! At least ',\
   max(param_col)+1,'columns needed for',mode)
   quit()

 n_f =[]
 for i in range(ff_nm_len):
  n_f.append(srl[i])
 ff_nm.append(n_f)

 n_p=[]
 for i,j in zip(param_col, param_typ):
  if j==0:
   k = int(srl[i])
  elif j==1:
   k = float(srl[i])
  n_p.append(k)
 ff_param.append(n_p)


def ChargeGaussian(fnm): #extract guassian esp charge from output
 c=[]
 f=open(fnm,'r')
 found=0
 for rl in f:
  if 'Charges from ESP fit' in rl:
   found=1
   break
 if found==0:
  print( 'Wrong guassian file')
  sys.exit()
 print( 'reading charges...')
 #print rl
 for rl in f:
  #print rl,
  if 'Charges from' in  rl:
   break
  srl=rl[:-1].split()
  if len(srl)==3:
   c.append(float(srl[2]))

 f.close()
 return c

#input nbfix parameters, LJ paramters [C12, C6]
#output list of atom type that has A A def in nbfixlist
#also Emin and Rmin as a  list 
def nbfix_self(nblist, nblj):
 slist=[]
 plist=[]
 for i,j in zip(nblist, nblj):
  if i[0]==i[1]:
   #print i[0], i[1], j[0], j[1]
   if i[0][0]=='H': continue
   c12=j[0]
   c6=j[1]
   epsilon = c6*c6/c12/4.
   sigma = (c12/c6)**(1/6.)
   #rmin=sigma
   rmin = sigma*(2.**(1/6.)) 
   if i[0] in slist:
    print( 'duplicate atom type', i[0])
    continue
   slist.append(i[0])
   plist.append([epsilon, rmin/2.])
 return (slist, plist)


#assign parameters to new types according to their analogy to mapped types
#for bond and angle, r0 or q0 will beused from pdb(heavy atoms) or qm(with H) geometry
#for improper and dihe, the existing values will be transfered
#FF will get the new parameters and will be changed
def assign_param_newtyp(FF, pdb, qm,conn_qm,typ, newtypmapping, ffatyp, Amber=True):
 conn=conn_qm
 conn.gen_impr_list(qm.atomnms,Amber)
 print( conn.impr)
 typnm=[i[0] for i in ffatyp]
 nbnm=[i[0] for i in FF.nbs]
 print( 'New atom types have been found')
#tranfer new nbond parameters
 for i in newtypmapping:
  print( i[0],'mapped to',i[1])
  FF.nbs.append([i[0]])
  idx=nbnm.index(i[1])
  
  FF.nb_params.append(FF.nb_params[idx])

#check if all the heavy atoms of qm also present in pdb
 nm_qm=qm.atomnms
 nm_pdb=pdb.atomnms
 xp_qm= qm.xp
 xp_pdb = pdb.xp
 for i in nm_qm:
  if i[0]!='H':
   if not i in nm_pdb:
    print( 'Err! ',i,'in QM geometry not in PDB geometry')
    sys.exit()

 #create a decoy of original type
 org_typ=[]
 maptyp=dict([(i[0],i[1]) for i in newtypmapping])
 for j in typ:
  if j in maptyp:
   org_typ.append(maptyp[j])
  else:
   org_typ.append(j)

 
#select all the entries of index multiplets involving with new atomtypes

 impr= pick_FF_newtyp(conn.impr, typ, maptyp)
 c12= pick_FF_newtyp(conn.c12, typ, maptyp)
 c13= pick_FF_newtyp(conn.c13, typ, maptyp)
 c14= pick_FF_newtyp(conn.c14, typ, maptyp)

# print impr

 print( 'Trying to set parameters for new types....')
#this function will change conn.impr
#so the conn_qm.impr will be changed after calling assign_param_newtyp
#only parameters involving new types will be scanned
#using org_typ, we first map new_types to old types
 mis=1

 while mis>0:
  cim_set2FF=[ match_impr_amber_list(i, FF.imprs, org_typ)   for i in impr]

 #sys.exit()

#the types mapped to original ones
  c12_typ_org = [get_nb_nm_list(i, c12, org_typ) for i in range(len(c12))]
  c13_typ_org = [get_nb_nm_list(i, c13, org_typ) for i in range(len(c13))]
  c14_typ_org = [get_nb_nm_list(i, c14, org_typ) for i in range(len(c14))]
  cim_typ_org = [get_nb_nm_list(i, impr, org_typ) for i in range(len(impr))]
#the new types
  c12_typ = [get_nb_nm_list(i, c12, typ) for i in range(len(c12))]
  c13_typ = [get_nb_nm_list(i, c13, typ) for i in range(len(c13))]
  c14_typ = [get_nb_nm_list(i, c14, typ) for i in range(len(c14))]
  cim_typ = [get_nb_nm_list(i, impr, typ) for i in range(len(impr))]

 
#nm used for new type
  c12_nm = [get_nb_nm_list(i, c12, qm.atomnms) for i in range(len(c12))]
  c13_nm = [get_nb_nm_list(i, c13, qm.atomnms) for i in range(len(c13))]
  c14_nm = [get_nb_nm_list(i, c14, qm.atomnms) for i in range(len(c14))]
  cim_nm = [get_nb_nm_list(i, impr, qm.atomnms) for i in range(len(impr))]

  c12_set2FF=[ match_FF_list(i, FF.bonds)   for i in c12_typ_org]
  c13_set2FF=[ match_FF_list(i, FF.angles)   for i in c13_typ_org]
  c14_set2FF=[ match_FF_list(i, FF.dihes)   for i in c14_typ_org]

# print FF.bonds
# print c12
# print c12_typ_org
# print c12_set2FF
# sys.exit()
#now cx_set2FF show the index of FF entry available for cx with mapping
#need to add new entry according to cx_set2FF
#for bond and angle, change r0 or q0, keep force constant

  nentry_bond=0
  nentry_angle=0
  nentry_dihe=0
  nentry_impr=0
  for inm, ityp, i2FF in zip(c12_nm, c12_typ, c12_set2FF):
   if (not ityp in FF.bonds) and len(i2FF)>0:
    #print ityp, FF.bonds[i2FF[0]]
    FF.bonds.append([i for i in ityp])
    FF.bond_params.append([i for i in FF.bond_params[i2FF[0]]])
    #correct r0
    if ityp[0][0]!='H' and ityp[1][0]!='H':
    #all heavy
     xp1 = xp_pdb[nm_pdb.index(inm[0])]
     xp2 = xp_pdb[nm_pdb.index(inm[1])]
    else:
     xp1 = xp_qm[nm_qm.index(inm[0])]
     xp2 = xp_qm[nm_qm.index(inm[1])]
    FF.bond_params[-1][1]= geom.bond(xp1, xp2)
    nentry_bond+=1
  for inm, ityp, i2FF in zip(c13_nm, c13_typ, c13_set2FF):
   if (not ityp in FF.angles) and len(i2FF)>0:
    FF.angles.append([i for i in ityp])
    FF.angle_params.append([i for i in FF.angle_params[i2FF[0]]])
    #correct q0
    if ityp[0][0]!='H' and ityp[1][0]!='H' and ityp[2][0]!='H':
    #all heavy
     xp1 = xp_pdb[nm_pdb.index(inm[0])]
     xp2 = xp_pdb[nm_pdb.index(inm[1])]
     xp3 = xp_pdb[nm_pdb.index(inm[2])]
    else:
     xp1 = xp_qm[nm_qm.index(inm[0])]
     xp2 = xp_qm[nm_qm.index(inm[1])]
     xp3 = xp_qm[nm_qm.index(inm[2])]
    FF.angle_params[-1][1]= geom.angle(xp1, xp2, xp3)
    nentry_angle+=1
  
#for dihe and improper, copy everything

  for inm, ityp, i2FF in zip(cim_nm, cim_typ, cim_set2FF):
   if (not ityp in FF.imprs) and len(i2FF)>0:
    FF.imprs.append([i for i in ityp])
    FF.impr_params.append([i for i in FF.impr_params[i2FF[0]]])
    nentry_impr+=1
#for dihe notice issue of multiplicity
  for inm, ityp, i2FF in zip(c14_nm, c14_typ, c14_set2FF):
   for j in i2FF:
    if is_dihe_in_FF(ityp, FF.dihe_params[j][1], FF.dihes, FF.dihe_params):
     continue
   
    FF.dihes.append([i for i in ityp])
    FF.dihe_params.append([i for i in FF.dihe_params[j]])
    nentry_dihe+=1

  print( 'New parameters set: BOND', nentry_bond,'ANGLE', nentry_angle,'IMPROPER',\
  nentry_impr,'DIHEDRAL',nentry_dihe)




  mis=0
  mis+=helper.show_missing(c12_set2FF, c12_typ_org, c12_nm,'BOND',ffatyp)
  mis+=helper.show_missing(c13_set2FF, c13_typ_org, c13_nm,'ANGLE',ffatyp)
  mis+=helper.show_missing(c14_set2FF, c14_typ_org, c14_nm,'DIHEDRAL',ffatyp)
  mis+=helper.show_missing(cim_set2FF, cim_typ_org, cim_nm,'IMPROPER',ffatyp)
  if mis>0:
   helper.input_param(FF, ffatyp)
#  break


#create new list of multiplet for any involvng types in select typ
#sl[][] and reutrn nl[][] nl[] point to sl[], 
#so if nl[] changes, sl[] also changes
def pick_FF_newtyp(sl, typ, seltyp):
 nl=[]
 for i in sl:
  for j in i:
   if typ[j] in seltyp:
    nl.append(i)
    break
 return nl 
   


