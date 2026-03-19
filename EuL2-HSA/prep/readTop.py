import top_db
import geom
import PDB
import connectivity
import helper
import sys
import TwoD_data
#from connectivity import *

def cap_str(s):
 re=''
 for i in s:
  if ord(i)>=ord('a') and ord(i)<=ord('z'):
   re+=chr(ord(i)-ord('a')+ord('A'))
  else:
   re+=i
 return re


def show_ic(IC):
 for i in IC:
  print( 'IC',end='')
  for j in range(0,4): print( '%5s'%(i[j]),end='')
  for j in range(4,9): print( '%10.4f'%(i[j]),end='')
  print()

def show_simple_list(s_list):
 for i in s_list: print (i,end='')
 print()

def show_zip_list(a,b):
 for i,j in zip(a,b):
  print( i,j)

#return any word in word_list that is the textline
#otherwise return None
def word_in_line(rl, w_list):
 re=None
 for w in w_list:
  if w in rl:
   re=w
   break
 return re 

#to see if any of a list of words
#could be found in the remaining of the text
#return the word that is found, otherwise return None
def find_next(f,w_list):
 re=None
 for rl in f:
  re=None
  for w in w_list:
   if w in rl:
    re=w
    break
  if not re is None:
   break
 #print re
 return re

def change_name(nms, nm_map):
 for i in range(len(nms)):
  old_nm = nms[i]
  nms[i]= nm_map[old_nm]

def change_name_list(nm_lists, nm_map):
 for i in nm_lists:
  change_name(i, nm_map)

#convert name atom to uniq 4-char name for PDB
#0.......9 A....Z
def convertName(nm_list):
 re=[]
 for n in nm_list:
  if len(n)>4:
   l=n[-1:]
   shi = ord(l)-ord('a')
   nn=n[0]
   for ki in range(1, 4):
    k=n[ki]
    if ord(k) >= ord('0') and ord(k)<ord('9'):
     cod = ord(k)-ord('0')
    else:
     cod = ord(k)-ord('A')+10
    if k==1: cod-=shi
    nn+= cod2char(cod)
   if nn in nm_list:
    print( 'Err!',nn,' is already in provided name list')
    quit()
  else:
   nn=n
  re.append(nn)  
 return re

def cod2char(cod):
 if cod<0 or cod>=36: 
  print( 'Err code')
  quit()
 a_cod = 36 - cod -1
 if a_cod<10:
  return chr(a_cod+ ord('0'))
 else:
  return chr(a_cod-10 + ord('A')) 

#deep copy list, two levels
def copy_list(lB):
 lA=[]
 for i in lB:
  temp=[]
  for j in i:
   temp.append(j)
  lA.append(temp)
 return lA

#deep copy dict, with value as one level list
def copy_dict(dB):
 dA={}
 for i in dB:
  dA[i]=dB[i][:]
#  print i, dA[i]
# print dA
 return dA

#delete a list of items where some item contains nm
def del_item(arg, nm):
 temp=[]
 for i in arg:
  if not nm in i:
   temp.append(i)
  else:
   print( i,'deleted' )
 del arg
 return temp
#merge Blist to Alist, name of elements in B will converted 
#according to dict t
def merge_list(lA, lB, t):
 for i in lB:
  temp=[]
  for j in i:
   if j in t:
    temp.append(t[j])
   else:
    temp.append(j)
  lA.append(temp)

def show_list(l,n):
 for i in l:
  print( n,':',i)

def showf_list(l,n,nw):
 for i,j in enumerate(l):
  if i%nw==0:
   print( n,'   ',end='')
  for k in j:
   print( k,'   ',end='')
  if i%nw==nw-1:
   print()
 print()

def show_dict(l,n):
 for i in l:
  print( n,':',i,'::',l[i])

def showf_dict(l,t,n,ncol):
# print l
# print t
 for i in l:
  print( n,'   ', i,'   ',end='')
  for kk,j in enumerate(t[i]):
   if kk<ncol:
    print( j,'   ',end='')
  print()


#residue topology class
class ResTop:
 """
 Data Struct for Res Top
 atoms_attr: dictionary of attrs of atoms with atom names as keys
 bonds: doublet for bonds
 imprs: quatrplet for impropers
 nextres: resname for next residue in the file
 """
#delete atom
 def del_atom(self, anm):
  if anm in self.atom_attr:
   print( anm,'deleted')
   self.bonds=del_item(self.bonds,anm)
   self.imprs=del_item(self.imprs,anm)
   del self.atom_attr[anm]
   self.atomnms.remove(anm)
  
 def del_atoms(self,al):
  for i in al:
   self.del_atom(i)
#check if atom exists
 def check_atomnm(self, arg):
  for i in arg:
   if i[0]=='-' or i[0]=='+':
    continue
   if not i in self.atomnms:
    print( i,'not found in known atoms')
    quit() 
#load content from list of string
 def load_struct(self,struct, arg, narg):
  if len(arg)==0:
   return
  if len(arg)%narg!=0:
   print( 'Following structural info is not correct:')
   print( arg)
   quit()
  temp=[]
  for i,j in enumerate(arg):
   temp.append(j)
   if len(temp)==narg:
    if temp[narg-1]<temp[0]:
     temp.reverse()
    if False:  #temp in struct:
     print( 'Warning!',temp,'already exist.')
    else:
     self.check_atomnm(temp)
     struct.append(temp)
    temp=[]
    
#load atom info
 def load_atom_gmx(self, srl, istop=False):
  n_srl=[]
  n_srl.append('')
  ss=srl[4]
  if istop:
   if not ss in self.__newnm:
    self.__newnm.append(ss)
    self.__newidx.append(0)
   idx = self.__newnm.index(ss)
   self.__newidx[idx]+=1
   ss+='_'+str(self.__newidx[idx])
  n_srl.append(ss)
  n_srl.append(srl[1])
  n_srl.append(srl[6])
  n_srl.append(srl[7])
  #print istop
  self.load_atom(n_srl, False, istop)

 def load_atom(self, srl, isCharmm=False, istop=False):
  if not isCharmm:
   anm, tnm, charge, mass= (srl[1], (srl[2]), float(srl[3]), float(srl[4]))
  else:
   #print srl
   anm, tnm, charge = (srl[1], cap_str(srl[2]), float(srl[3]))
   mass=0.0
  #print istop
  if anm in self.atomnms and (not istop):
   print( 'Duplicated atom names',anm)
   quit()
  self.atoms.append([tnm, charge, mass])
  self.atomnms.append(anm)
#load bond info
 def load_bond(self, srl):
  arg=srl[1:]
  if len(arg)>0:
   self.load_struct(self.bonds, arg, 2)
#load improper info
 def load_impr(self, srl):
  arg=srl[1:]
  if len(arg)>0:
   self.load_struct(self.imprs, arg, 4)

#read a line in GMX format
#which way the data is loaded is according to term in load_rule
 def read_gmx_data(self, srl, term, load_rule, exemptlist=[]):
  if term in exemptlist: return
  #print term, srl
  #print load_rule
  
  l_rule = load_rule[term]
  l_type_pos = l_rule[4]
  l_subtype= l_rule[5]
  ta = int(l_type_pos)
  #if l_type_pos!=None:
  #at least type column should be present in data line
  if ta >= helper.len_comment(srl):
   print( srl,'has not enough terms for',term)
   sys.exit()
  tb = int(srl[ta])
  if tb in l_subtype:
   self.read_gmx_data(srl, l_subtype[tb] ,load_rule)
   return
  #load params field
  param_obj = l_rule[2]
  param_col_field= l_rule[3]
  max_col = max(param_col_field)
 # if max_col>=len(srl):
 #  print srl,'has not enough terms for',term
 #  sys.exit()
  t_srl=[]
  #if specified params are persent, read them, otherwise append []
  if max_col<helper.len_comment(srl):
   for i in param_col_field:
    t_srl.append(srl[i])
  param_obj.append(t_srl)

  #load name field
  name_obj = l_rule[0]
  name_field_col_idx =l_rule[1]
  n_srl=[]
  natom = len(self.atomnms)
  for i in name_field_col_idx:
   ok=True
   if not PDB.is_idx(srl[i]): ok=False
   idx=int(srl[i])
   if idx<1 or idx>natom: ok=False
   if not ok:
    print( 'Err in reading',srl,'for',term)
    quit()
   n_srl.append(self.atomnms[idx-1])
  self.load_struct(name_obj, n_srl, len(name_field_col_idx))

#read topology from GMX file
 def readGMX(self,f, opt,readCon=False, istop=False):
  if not readCon:
   lf =find_next(f,['moleculetype'])
   if lf is None:
    print( 'Err! Wrong GMX topology format.')
    quit()
  self.resnm = None
  self.nextres = False

  while self.resnm is None:
   
   for rl in f:
    srl=rl[:-1].split()
    if len(srl)==0: continue
    if srl[0][0]==';' : continue
    #print srl
    if opt=='ALL' or srl[0]==opt:
     self.resnm=srl[0]
     break
   #print self.resnm,'#' 
   if  self.resnm is None:
    lf =find_next(f,['moleculetype'])
    if lf is None:
     print( 'Err! Wrong GMX topology format.')
     quit()


  lf = find_next(f,['atoms'])
  if lf is None:
   print( 'Err! Cannot find atom section')
   quit()
  self.atomnms=[]
  self.atoms=[]
  self.bonds=[]
  self.bond_params=[]
  self.angles=[]
  self.angle_params=[]
  self.dihes=[]
  self.dihe_params=[]
  self.imprs=[]
  self.impr_params=[]
  self.pairs=[]
  self.pair_params=[]
  self.exclusions=[]
  self.exclusion_params=[]
  self.IC=[]
  #load atoms
  #                     field for names         field for params
  load_rule = top_db.top_gmx_load_rules(self)
  terms= load_rule.keys()
  terms.append('moleculetype')
  w=None
  
  if istop:
   self.__newnm=[]
   self.__newidx=[]
  for rl in f:
   srl=rl[:-1].split()
   if len(srl)==0: continue
   if srl[0][0]==';': continue
   w = word_in_line(rl, terms)
   if not w is None:
    break
   self.load_atom_gmx(srl,istop)

  if self.atoms==[]:
   print( 'Err! No atom found')
   quit() 
  #if istop:
  # self.atom_attr = self.atoms
  #else:   
  self.atom_attr = dict([(x,y) for x,y in zip(self.atomnms, self.atoms)])  
  
  if w is None:
  #this means file ends
   return
  #readMode=terms.index(w)
  #if readMode==len(terms)-1: return
  if w=='moleculetype':
  #there are more residues after this one
   self.nextres = True 
   return
  for rl in f:
   #print w, rl,
   srl=rl[:-1].split()
   if len(srl)==0: continue
   if srl[0][0] in [';','#','!']: continue
   w_t = word_in_line(rl, terms+['position_restraints','system','molecules'])
   if not w_t is None:
    #readMode=terms.index(w)
    #if readMode==len(terms)-1: return
    if w_t=='moleculetype':
     self.nextres=True 
     return
    w= w_t
    continue
   self.read_gmx_data(srl, w, load_rule,\
   ['position_restraints','system','molecules']) #need to continue here

#constructor from rtf file, f is a file ref
#istop=True to load a complete topology file for GMX
 def __init__(self, f=None, rtB=None, opt='ALL', readCon=False,\
              readGMX=False,istop=False):
  if readGMX:
   self.readGMX(f,opt,  readCon, istop)
   return
  if not rtB==None:
   self.copy(rtB)
   return
  found = False
  self.resnm = opt
  self.atomnms=[]
  if not readCon:
   self.atoms=[]
   self.bonds=[]
   self.bond_params=[]
   self.angles=[]
   self.angle_params=[]
   self.dihes=[]
   self.dihe_params=[]
   self.pairs=[]
   self.pair_params=[]
   self.imprs=[]
   self.impr_params=[]
   self.IC=[]
   if f==None: return
   for rl in f:
    if ('RESI' in rl or 'PRES' in rl)and (opt=='ALL' or opt in rl):
     found = True
     srl=rl[:-1].split()
     self.resnm=srl[1]
     if 'PRES' in rl: self.resnm+='_PRES'
     break
   if not found:
    print( 'Err! Residue',opt,'not found.')
    quit()
  else:
   self.resnm=opt
  self.atoms=[]
  self.bonds=[]
  self.bond_params=[]
  self.angles=[]
  self.angle_params=[]
  self.dihes=[]
  self.dihe_params=[]
  self.imprs=[]
  self.pairs=[]
  self.pair_params=[]
  self.IC=[]
  self.impr_params=[]
  self.nextres = 'NONE'
  if f==None: return
  #if 'RESI' in rl:
  #self.resnm = rl[:-1].split()srl[1]
  for rl in f:
#   print rl,
   srl=rl[:-1].split()
   if not '_PRES' in self.resnm:
    if 'DELETE' in rl: continue
    if 'ATOM' in rl:
     self.load_atom(srl, True)
    if 'BOND' in rl or 'DOUBLE' in rl:
     self.load_bond(srl)
    if 'IMPR' in rl:
     self.load_impr(srl)
#   if 'PATCH' in rl:
#    break
   if 'RESI' in rl or 'PRES' in rl:
    self.nextres = srl[1]
    if 'PRES' in rl: self.nextres+='_PRES'
    break
  self.atom_attr = dict([(x,y) for x,y in zip(self.atomnms, self.atoms)])
#warning! don't merge with self
 def merge(self, rtB, suf='a'):
  t={}
  for i in rtB.atomnms:
   tnm=i
   if i in self.atomnms:
    tnm+=suf
   t[i]=tnm
  for i in rtB.atomnms:
   self.atomnms.append(t[i])
   self.atom_attr[t[i]]=rtB.atom_attr[i][:]
  merge_list(self.bonds, rtB.bonds, t)
  merge_list(self.imprs, rtB.imprs, t)  

 def tcharge(self):
  t=0.0;
  for i in self.atom_attr:
   t+=self.atom_attr[i][1]
  return t

 def gettypetuple(self, flist):
  re=[]
  for i in flist:
   #print self.atom_attr
   if not i in self.atom_attr: print( i,'not found in top')
   re.append(self.atom_attr[i][0])
  return re

 def showf(self):
  show_list(self.atomnms,'NAME')
  show_dict(self.atom_attr,'ATTR')
  show_list(self.bonds, 'BOND')
  show_list(self.imprs,'IMPR')   
  
 def show(self):
  print( 'RESI   ',self.resnm,'   {0:.2f}'.format(self.tcharge()))
  showf_dict(self.atomnms, self.atom_attr,'ATOM',2)
  showf_list(TwoD_data.unq_set(self.bonds,False),'BOND',4)
  showf_list(TwoD_data.unq_set(self.imprs,False),'IMPR',2)
  showf_list(TwoD_data.unq_set(self.angles,False),'ANGLE',3)
  showf_list(TwoD_data.unq_set(self.dihes,False),'DIHE',2)
  if not self.IC is None:
   show_ic(self.IC)
  print( 'PATCHING FIRST NONE LAST NONE')


 def show_param(self):
  show_zip_list(self.bonds,self.bond_params)
  show_zip_list(self.angles,self.angle_params)
  show_zip_list(self.imprs,self.impr_params)
  show_zip_list(self.dihes,self.dihe_params)

#deep copy
 def copy(self,rtB):
  self.resnm = rtB.resnm
  self.nextres = rtB.nextres
  self.bonds=copy_list(rtB.bonds)
  self.imprs=copy_list(rtB.imprs)
  self.atomnms=rtB.atomnms[:]
  self.atom_attr=copy_dict(rtB.atom_attr)

#make all atom name not more than 4 char
 def rename(self):
  s_nms = convertName(self.atomnms)
  nm_map = dict([(i,j) for i,j in zip(self.atomnms, s_nms)])
  #print nm_map
  n_attr= dict([(nm_map[i], self.atom_attr[i]) for i in self.atomnms])
  change_name(self.atomnms, nm_map)
  change_name_list(self.bonds, nm_map)
  change_name_list(self.imprs, nm_map)
  self.atom_attr = n_attr

 def add_atom(self,anm, tnm, c):
  if anm in self.atomnms:
   print( anm,'already exist')
   return
  self.atomnms.append(anm)
  self.atom_attr[anm]=[tnm,c]

 def add_bond(self,anm,bnm):
  if not anm in self.atomnms or not bnm in self.atomnms:
   print( anm,'or',bnm,'doesnot exist')
   return
  anm,bnm=(min(anm,bnm),max(anm,bnm))
  if [anm,bnm] in self.bonds:
   print( 'Bond',anm,'-',bnm,'exists')
   return
  self.bonds.append([anm,bnm])

 def integrity_check(self):
  construct_conn(self,'top')
  construct_atomSig(self,None,False)
  if len(self.atomSigs[0])<len(self.atomnms)*3:
   print( 'Warning! Some part of molecule is disconnected from the other.')

 def update_conn(self, maskoff=None):
  construct_conn(self,'top')
  if maskoff==None:
   construct_atomSig(self)
  else:
   maskofflist=[]
   for i in maskoff:
    if i in self.atomnms:
     maskofflist.append(self.atomnms.index(i))
   construct_atomSig(self, maskofflist)

  sort_sig(self.atomSigs)
  construct_uniq(self)

 #make internal coordinate according to xp list
 def make_IC(self, xp):
  connTop= connectivity.connect(None, True,self)
  #print len(connTop.c12)
  #print len(connTop.c13)
  #print len(connTop.c14)
  conn= connTop.get_direct_conn()
  self.IC =[]
  
  n=len(self.atomnms)
  atomnms=self.atomnms
  mark=[0 for i in range(n)]
  #find heavy atom has max neibor connectivity to other heavy atoms
  maxcon=[]
  for ai in range(n):
   m=0
   for i in  conn[ai]:
   
    if atomnms[i][0]!='H':
     m+=1
     for j in conn[i]:
      if j!=ai and atomnms[j][0]!='H':
       m+=1
   maxcon.append(m)
  if max(maxcon)<2 : return
  first_try=maxcon.index(max(maxcon))
  t_l=[]
  for i in conn[first_try]:
   m=0
   for j in conn[i]:
    if atomnms[j]!='H': m+=1
   t_l.append(m)

  idx = t_l.index(max(t_l))
  t_l[idx]=-1
  sec_try=conn[first_try][idx]

  q=[first_try, sec_try]
  q_hist=[[sec_try, first_try],[first_try, sec_try]]
  mark[first_try]=1
  mark[sec_try]=1
  head=0
  while head<len(q):

   c_i = q[head]
   h_i=q_hist[head]
   #print c_i,atomnms[c_i], h_i
   #at current point, only one of the connected points
   #will be represented by dihe, other by impr
   #this one should be max connection to other heavy atom
   t_l=[]
   for i in conn[c_i]: 
     
     m=0
     if atomnms[i][0]!='H': m+=1
     for j in conn[i]:
      if atomnms[j][0]!='H': m+=1
     if mark[i]==1: m=-1
     t_l.append(m)
   #if it is H, m=0
   #if it has no unmarked atoms, m all -1

   #if no further connection
   #print conn[c_i], t_l,
   #for k in conn[c_i]: print atomnms[k],mark[k],
   #print

   if len(t_l)==0:
    head+=1
    print( 'back')
    continue
   #if certain points have been marked
   #it can be used as for impr construction
   idx_impr=-1
   last_visit=h_i[-2]
   idx_last=conn[c_i].index(last_visit)
   t_l[idx_last]=100 
   if min(t_l)<0:
    idx_impr=conn[c_i][t_l.index(min(t_l))]
   t_l[idx_last]=-1
   idx_further=-1    
   #if there is one point to go further
   if max(t_l)>=0:
    idx_further = conn[c_i][t_l.index(max(t_l))]
    #mark[idx_further]=1
    
    if idx_impr<0:
     #it used for other impr ic
     #itself used with dihe ic 
     idx_impr=idx_further
     ic_by_dihe(h_i+[idx_further],xp, atomnms,self.IC)
    else:
     #itself used with impr ic
     ic_by_impr(h_i+[idx_impr, idx_further], xp, atomnms,self.IC)
   #go over other unmarked and not to go further
   for i in conn[c_i]:
    if mark[i]==0 and i!=idx_further:
     ic_by_impr(h_i+[idx_impr, i], xp, atomnms,self.IC)
     
   #go further
   for i in conn[c_i]:
    if mark[i]==0:
     #print i, atomnms[i]
     mark[i]=1
     ss=h_i+[i]
     q.append(i)
     q_hist.append(ss)
   head+=1
    


  #then crowlin the whole molecules


#try to assign atom in top to atom in PDB

#construct IC by dihe def
def ic_by_dihe(l,xp, atomnms,IC):
 if len(l)<4: return
 nl=l[-4:]
 t=[atomnms[nl[0]],\
    atomnms[nl[1]],\
    atomnms[nl[2]],\
    atomnms[nl[3]],\
    geom.bond(xp[nl[0]],xp[nl[1]]),\
    geom.angle(xp[nl[0]],xp[nl[1]],xp[nl[2]]),\
    geom.dihedral(xp[nl[0]],xp[nl[1]],xp[nl[2]],xp[nl[3]]),\
    geom.angle(xp[nl[1]],xp[nl[2]],xp[nl[3]]),\
    geom.bond(xp[nl[2]],xp[nl[3]])]
 IC.append(t)

#construct IC by impr def
def ic_by_impr(l,xp, atomnms,IC):
 if len(l)<4: return
 nl=l[-4:]
 t=[atomnms[nl[2]],\
    atomnms[nl[0]],\
    '*'+atomnms[nl[1]],\
    atomnms[nl[3]],\
    geom.bond(xp[nl[2]],xp[nl[1]]),\
    geom.angle(xp[nl[2]],xp[nl[1]],xp[nl[0]]),\
    geom.dihedral(xp[nl[2]],xp[nl[0]],xp[nl[1]],xp[nl[3]]),\
    geom.angle(xp[nl[0]],xp[nl[1]],xp[nl[3]]),\
    geom.bond(xp[nl[1]],xp[nl[3]])]
 IC.append(t)



#================
#Construct uniqoueness of list
#add a list for each atom: unq[[idx of first ele that is the same as current, number of the same ele]]
#===============

def construct_uniq(obj):
 Sigs=obj.atomSigs
 unq=[]
 ref=0
 nrep=0
 for i,j in enumerate(Sigs):
  if sigCmp(j, Sigs[ref])==0:
   nrep+=1
   unq.append([ref,0])
  else:
   unq[ref][1]=nrep
   nrep=1
   ref=i
   unq.append([ref,0])

 unq[ref][1]=nrep
 obj.unq=unq
   
#===================================
#To construct geometry
#according to cnnoectivity
#adding the connective properties to given objs
#specifically, it will add two attributes to obj
#anm_ini[] , list of initals of atom names
#conn[[]], the list of indexs of connecting atoms for each atom
#in ascending order
#===================================
def construct_conn(obj, source,atoms= None ):
 if (source=='pdb'):
  construct_conn_pdb(atoms, obj)
 elif (source=='top'):
  construct_conn_top(obj)
 else: 
  print( 'Err! No such option %s for construction of connectivity'%source)
  quit()

def init_anm_conn(atoms):
 anm_ini=[]
#atom name initials
 for  i in atoms:
  anm_ini.append(i[1][0])
 conn=[]
 for i in range(len(atoms)):
  conn.append([])
 return anm_ini, conn

def construct_conn_pdb(atoms, obj):
 anm_ini, conn = init_anm_conn(atoms)
 n = len (atoms)
 for i in range(n-1):
  atom_xi = PDB.load_pdb_xyz(atoms[i])
  nm_i = anm_ini[i] 
  for j in range(i+1,n):
   atom_xj = PDB.load_pdb_xyz(atoms[j])
   nm_j = anm_ini[j]
   if nm_i=='H' or nm_j=='H':
    cut = 1.3
   elif nm_i=='S' and nm_j=='S':
    cut=2.2
   else:
    cut = 1.9
   if geom.dist(atom_xi, atom_xj)<cut:
    conn[i].append(j)
    conn[j].append(i)
#both conn[i] and [j] should be sorted by this way
 obj.conn = conn
 obj.anm_ini = anm_ini

def construct_conn_top(obj):
 
 n = len (obj.atomnms)
 anm_ini=[]
 conn=[]
 for i in obj.atomnms:
  anm_ini.append(i[0])
  conn.append([])

 for ib in obj.bonds:
   anm=ib[0]
   bnm=ib[1]
   if not anm in obj.atomnms or not bnm in obj.atomnms:
    print( 'Err! some atom in bonds are not listed in atom name list.')
    quit()
   i = obj.atomnms.index(anm)
   j = obj.atomnms.index(bnm)
   nm_i = anm[0]
   nm_j = bnm[0]
   conn[i].append(j)
   conn[j].append(i)   
#conn[i] and [j] may not be sorted
 for i in conn:
  i.sort()

 obj.conn = conn
 obj.anm_ini = anm_ini


#================================
#construct atom signature
#which is a list [ini_connecting atom, num for generation,....] for each atom
#require connectivity (conn) and name initials (anm_ini) from obj
#defalut expoloration depth is 4
#================================

def construct_atomSig(obj, maskoff=None,withoutH=True,gener=1000):
 n=len(obj.anm_ini)
 conn = obj.conn
 anm_ini=obj.anm_ini
 #print conn
 atomSigs=[]
 for i in range(n):
#construct atomSig atom by atom
  if withoutH and anm_ini[i]=='H': continue
  if maskoff!=None:
   if i in maskoff: continue
  Sig=[]
  add_sig (Sig, i, conn, anm_ini,  maskoff, gener, withoutH)
  atomSigs.append(Sig)
  del Sig
 obj.atomSigs = atomSigs

def add_sig(Sig,i, conn, anm_ini, maskoff, gener, withoutH):
 #print i
 n=len(anm_ini)
 mark=[0 for j in range(n)]
 que=[]
 que_level=[]
 que.append(i)
 que_level.append(0)
 head=0
 mark[i]=1
#mask off designated atoms
 if maskoff!=None:
  for mi in maskoff: mark[mi]=1

 #print i,que
 while head<len(que):
  i_c= que[head]
  i_lev=que_level[head]
  #print head,':',i_c,'|',
  Sig.append(anm_ini[i_c])
  Sig.append(i_lev)
  Sig.append(i_c)
  if i_lev>=gener:
   head+=1
   continue
  #print i_c,':',conn[i_c]
  for k in conn[i_c]:
   #print k,
   if mark[k]==0 and not (anm_ini[k]=='H' and withoutH):
    mark[k]=1
    que.append(k)
    que_level.append(i_lev+1)
  head+=1
 

#================
#sort atom signiture
#non-recursive mergesort
#===========

def sigCmp(a,b):
 na=len(a)
 nb=len(b)
 if na<nb: return -1
 if na>nb: return 1
 for i in range(na):
  if i%3!=2:
   if a[i]<b[i]: return -1
   if a[i]>b[i]: return 1
 return 0


def sort_short(l, r, Sigs):
 if r-l==1:
  if sigCmp(Sigs[l], Sigs[r])>0: #sigCmp can be changed for sorting other data
   swp=Sigs[r]
   Sigs[r]=Sigs[l]
   Sigs[l]=swp

def merge_sort(l1,r1,l2,r2, Sigs):
 if r1!=l2-1:
  print( 'Two segments [',l1,r1,'][',l2,r2,'] are not continoues.')
  quit()
 La = [Sigs[i] for i in range(l1, r1+1)]
 Ra = [Sigs[i] for i in range(l2, r2+1)]
 nL = len(La)
 nR = len(Ra)
 s1=0
 s2=0
 for i in range(l1, r2+1):
   if s1<nL and s2<nR:
    if sigCmp(La[s1],Ra[s2])<0:
     Sigs[i]=La[s1]
     s1+=1
    else:
     Sigs[i]=Ra[s2]
     s2+=1
   elif s1<nL:
    Sigs[i]=La[s1]
    s1+=1
   else:
    Sigs[i]=Ra[s2]
    s2+=1
     

def sort_sig(Sigs):
 stk=[]
 l, r=(0, len(Sigs)-1)
 stk.append((l,r))
 head=0
#divid
 while (head<len(stk)):
  l, r = stk[head]
  if (r-l>1):
   m=int((r+l)/2)
   stk.append((l,m))
   stk.append((m+1,r))
  head+=1
 eo=1
# print stk
#conquer
 while (len(stk)>0):
  #print eo, stk
  if eo==1:
   le, re=stk.pop()
   sort_short(le, re, Sigs)
   eo=0
   continue
  if eo==0:
   lo, ro=stk.pop()
   sort_short(lo, ro, Sigs)
   merge_sort(lo, ro, le, re, Sigs)
   eo=1
   
def load_tops_gmx(fnm):
 tops=[]
 f=open(fnm,'r')
 n_res= ResTop(f, None, 'ALL', False, True)
 tops.append(n_res)
 while n_res.nextres:
  n_res=ResTop(f, None, 'ALL', True, True)
  tops.append(n_res)
 f.close() 
 return tops

def load_tops_rtf(fnm):
 tops=[]
 f=open(fnm,'r')
 n_res= ResTop(f, None, 'ALL', False, False)
 tops.append(n_res)
 #print n_res.resnm
 while n_res.nextres!='NONE':
  n_res=ResTop(f, None, n_res.nextres, True, False)
  tops.append(n_res)
  #print n_res.resnm
 f.close()
 return tops

def load_gmx_full_top(fnm):
 f=open(fnm,'r')
 top= ResTop(f, None, 'ALL', False,\
                    True, True)
 f.close()
 return top

def summary_atom_in_top(tops):
 re=[]
 mass=[]
 for itop in tops:
  anms=itop.atomnms
  attr=itop.atom_attr
  for i in anms:
   c_t = attr[i][0]
   if not c_t in re:
    re.append(c_t)
    mass.append(attr[i][2])
 return (re, mass) 

def load_atom_type(fnm, charmm=True):
 f=open(fnm,'r')
 t=[]
 
 for rl in f:
  if 'MASS' in rl:
   srl=rl[:-1].split()
   comm=''
   k=0
   for i in srl:
    if '!' in i:
     k=1
    if k==1: comm+=i+' ' 
   t.append([srl[2],float(srl[3]), comm])
 f.close()

 return t
