import sys
import readFF
import readTop
import geom
import numpy as np

substitute_head = ['bonds', 'exclusions', 'pairs', 'angles', 'dihedrals','cmap']
substitute_nitem = [2, 2, 2, 3, 4, 5]

guess_mass = {'C':12.0, 'H':5.0,'O':16.0,'N':14.0,'S':32.0}

aro_nm = ['CZ','CF','NZ']

class FFgmx:
# a simplified gmx FF class for output only
# only work with atomtype
 def __init__(self):
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
  self.gmx_atom_rmks=[]
  self.gmx_nb_rmks=[]
  self.gmx_pair_rmks = []
  self.gmx_bond_rmks = []
  self.gmx_angle_rmks = []
  self.gmx_dihe_rmks=[]

 def load_file(self,fnm): 
   self.txt = []
   self.sechead = []
   self.secpos = []
   self.warn = 0
   self.note = 0

   txt_string = [rl for rl in open(fnm,'r')]

   for _i,rl in enumerate(txt_string):
    self.txt.append(rl)
#    print(rl,end='')
    srl = rl[:-1].split()
    if len(srl)==0: continue
    if '[' in srl[0]:
     _rl1 = rl[:-1].replace('[',' ')
     _rl2 = _rl1.replace(']',' ')
     self.sechead.append(_rl2.split()[0])
     self.secpos.append(_i)

   self.secpos.append(len(txt_string))
   
# load atoms
   if 'atomtypes' in self.sechead:
    self.gmx_atom_rmks=[]
    _secidx = self.sechead.index('atomtypes')
    _txt = self.txt[self.secpos[_secidx]:self.secpos[_secidx+1]]
#    print(_txt)
    for _rl in _txt:
     srl = _rl[:-1].split()
     if _is_nondata_line(_rl): continue 
     self.gmx_atoms.append([srl[0]])
     self.gmx_atom_params.append([\
      float(srl[1]),float(srl[2]),float(srl[4]),float(srl[5])])
     self.gmx_atom_rmks.append('')

   if 'nonbond-params' in self.sechead:
    self.gmx_nb_rmks = []
    _secidx = self.sechead.index('nonbond-params')
    _txt = self.txt[self.secpos[_secidx]:self.secpos[_secidx+1]]
    for _rl in _txt:
     srl = _rl[:-1].split()
     if _is_nondata_line(_rl): continue
     self.gmx_nbs.append(sorted([srl[0], srl[1]]))
     self.gmx_nb_params.append([\
      float(srl[3]),float(srl[4])])
     self.gmx_nb_rmks.append('')


   if 'pairtypes' in self.sechead:
    self.gmx_pair_rmks = []
    _secidx = self.sechead.index('pairtypes')
    _txt = self.txt[self.secpos[_secidx]:self.secpos[_secidx+1]]
    for _rl in _txt:
     srl = _rl[:-1].split()
     if _is_nondata_line(_rl): continue
     self.gmx_pairs.append(sorted([srl[0], srl[1]]))
     self.gmx_pair_params.append([\
      float(srl[3]),float(srl[4])])
     self.gmx_pair_rmks.append('')

   if 'bondtypes' in self.sechead:
    self.gmx_bond_rmks = []
    _secidx = self.sechead.index('bondtypes')
    _txt = self.txt[self.secpos[_secidx]:self.secpos[_secidx+1]]
    for _rl in _txt:
     srl = _rl[:-1].split()
     if _is_nondata_line(_rl): continue
     self.gmx_bonds.append(sorted([srl[0], srl[1]]))
     self.gmx_bond_params.append([\
      float(srl[3]),float(srl[4])])
     self.gmx_bond_rmks.append('')

   if 'angletypes' in self.sechead:
    self.gmx_angle_rmks = []
    _secidx = self.sechead.index('angletypes')
    _txt = self.txt[self.secpos[_secidx]:self.secpos[_secidx+1]]
#    print('check angle') 
    for _rl in _txt:
#     print(_rl)
     srl = _rl[:-1].split()
     if _is_nondata_line(_rl): continue
     self.gmx_angles.append([srl[0], srl[1], srl[2]])
     self.gmx_angle_params.append([\
      float(srl[4]),float(srl[5])])
     self.gmx_angle_rmks.append('')


   if 'dihedraltypes' in self.sechead:
    self.gmx_dihe_rmks = []
    _secidx = self.sechead.index('dihedraltypes')
    _txt = self.txt[self.secpos[_secidx]:self.secpos[_secidx+1]]
    for _rl in _txt:
     srl = _rl[:-1].split()
     if _is_nondata_line(_rl): continue
     self.gmx_dihes.append([srl[0], srl[1]])
     self.gmx_dihe_params.append([\
      float(srl[3]),float(srl[4]), int(srl[5])])
     self.gmx_dihe_rmks.append('')


 def show(self,fnm):
  f = open(fnm,'w')

  f.write('[ atomtypes ]\n')
  for _i,_j,_k in zip(self.gmx_atoms,\
                   self.gmx_atom_params, self.gmx_atom_rmks):
    f.write('%s   %f   %f   A   %g   %g ; %s\n'\
     %(_i[0], _j[0],_j[1],_j[2],_j[3], _k))

  f.write('\n')
  f.write('[ nonbond-params ]\n')
  for _i,_j,_k in zip(self.gmx_nbs,\
                   self.gmx_nb_params, self.gmx_nb_rmks):
    f.write('%s   %s   1      %g   %g ; %s\n'\
     %(_i[0], _i[1],_j[0],_j[1], _k))
  f.write('\n')
  f.write('[ pairtypes ]\n')
#  print(self.gmx_pair_rmks)
  for _i,_j,_k in zip(self.gmx_pairs,\
                   self.gmx_pair_params, self.gmx_pair_rmks):
#    print(_i,_j,_k)
    f.write('%s   %s   1      %g   %g ; %s\n'\
     %(_i[0], _i[1],_j[0],_j[1], _k))

  f.write('\n')

  f.write('[ bondtypes ]\n')
  for _i,_j,_k in zip(self.gmx_bonds,\
                   self.gmx_bond_params, self.gmx_bond_rmks):
    f.write('%s   %s   1      %f   %f ; %s\n'\
     %(_i[0], _i[1],_j[0],_j[1], _k))

  f.write('\n')

#  print(self.gmx_angles)
  f.write('[ angletypes ]\n')
  for _i,_j,_k in zip(self.gmx_angles,\
                   self.gmx_angle_params, self.gmx_angle_rmks):
    f.write('%s   %s   %s  1      %f   %f ; %s\n'\
     %(_i[0], _i[1],_i[2],_j[0],_j[1], _k))

  f.write('\n')

  f.write('[ dihedraltypes ]\n')
  for _i,_j,_k in zip(self.gmx_dihes,\
                   self.gmx_dihe_params, self.gmx_dihe_rmks):
    f.write('%s   %s    1      %f   %f  %d ; %s \n'\
     %(_i[0], _i[1],_j[0],_j[1],_j[2], _k))

  f.write('\n')


  f.close() 

class mod_typ:

 def __init__(self, anm, atyp, atypB):
  self.anm=anm
  self.atyp = atyp
  self.atypB = atypB

 def lookup(self, _a):
  if not _a in self.anm:
   return (None, None)
  else:
   _id = self.anm.index(_a)
   return (self.atyp[_id], self.atypB[_id])

def _is_nondata_line(rl):
 srl = rl[:-1].split()
 if len(srl)==0: return True
 if srl[0][0] == ';' or srl[0][0] == '#' or srl[0][0] == '[': return True
 return False


def find_end_idx_sec(rl, nitem):
# find the starting position in a line for non-index txt

 _is_space=True

 _ic = 0
 _st_pos = -1
 for _i in range(len(rl)):
#  print(_ic,'|%s|'%rl[_i])
  if _is_space:
   if rl[_i]!=' ' and rl[_i]!='\t' and rl[_i]!='\n':
    _is_space = False
  elif rl[_i]==' ' or rl[_i]=='\t' or rl[_i]=='\n':
    _ic+=1
    _is_space = True
    if _ic== nitem:
     _st_pos = _i
     break

 if _st_pos<0:
  print('Error: the following line does not contain',nitem,'index items')
  print('Error:',rl,end='')
  sys.exit()
    
 return _st_pos

class topF:

 def __init__(self, fnm=None, txt_string=None):

   self.txt = []
   self.sechead = []
   self.secpos = []
   self.warn = 0
   self.note = 0

   if txt_string is None:
    txt_string = [rl for rl in open(fnm,'r')]

   for _i,rl in enumerate(txt_string):
    self.txt.append(rl)
#    print(rl,end='')
    srl = rl[:-1].split()
    if len(srl)==0: continue
    if '[' in srl[0]:
     _rl1 = rl[:-1].replace('[',' ')
     _rl2 = _rl1.replace(']',' ')
     self.sechead.append(_rl2.split()[0])
     self.secpos.append(_i)

   self.secpos.append(len(txt_string))
   
#   print(self.sechead)
   _sec_idx = self.sechead.index('atoms')
   _load_atom_line = self.txt[self.secpos[_sec_idx]: self.secpos[_sec_idx+1]]

   self.atoms_update(_load_atom_line)


 def atoms_update(self, _load_atom_line):

   self.aid = []
   self.atyp = []
   self.atomtxt=[]  
   self.atomnm=[]
   self.resnm = []
   self.resid =[]
   self.cgrp = []
   self.resid_unq = []
   _last_rid = -1

#   print(_load_atom_line)
   for rl in _load_atom_line:
    if _is_nondata_line(rl): continue

#    print(rl,end='')
    self.atomtxt.append(rl)
    srl=rl[:-1].split()
    self.atomnm.append(srl[4])
    self.resnm.append(srl[3])
    self.atyp.append(srl[1])
    self.cgrp.append(int(srl[5]))
    self.aid.append(int(srl[0]))
    self.resid.append(int(srl[2]))
    if _last_rid!=int(srl[2]):
     self.resid_unq.append(int(srl[2]))
     _last_rid = int(srl[2])

 def add_atoms(self, atoms_added,mod_check,_r_st,_r_end, FF,\
               debug=False):
 # atoms_added atoms to be added
 # _r_st, _r_end : starting and ending pos of mod residue
 # mod_check: atoms whose types may be changed in this residue
 # FF : force filed obj
 # this code only update atoms of topology by adding new atoms
 # it also does parameter check

  gmx_atom = [_i[0] for _i in FF.gmx_atoms]
  gmx_atom_p = FF.gmx_atom_params

  _nadd =len(atoms_added)-3 

# sanity check
# assuming old topoloy has its parameters well defined
  if debug:
   print(self.atyp)
#   print(atoms_added)
   print(mod_check.atyp)

  for _i in self.atyp:
   if not _i in gmx_atom:
    print('Error: atom type',_i,'not defined in force field.') 
    sys.exit()

# all the info for each atom stored in a text line     
  _new_txt=[]

  _natom = len(self.atomtxt)

  _to_insert = True

  for _i in range(_natom+1):

   if _i< _r_st:
    _new_txt.append(self.atomtxt[_i])
    continue

   if _i>= _r_end:

    if _to_insert:
# insert new section here
     _to_insert=False
     aid_b = self.aid[_i-1]
     resid = self.resid[_i-1]
     resnm = self.resnm[_i-1]
     cgrp = self.cgrp[_i-1]
     _cnt = 0
     for _atom_i in atoms_added[3:]:
      _anm = _atom_i[1]
      _atyp, _atypB = mod_check.lookup(_anm)
      _cnt+=1
      if _atyp in gmx_atom:
       _id = gmx_atom.index(_atyp)
       _par = gmx_atom_p[_id]
       _new_txt.append('%d  %s  %d  %s  %s  %d  %f  %f'%\
       (aid_b+_cnt,_atyp,\
        resid, resnm,\
        _anm, cgrp,\
        _par[1], _par[0])\
        +' ; added atom\n')
      elif _atypB in gmx_atom:
       _id = gmx_atom.index(_atypB)
       _par = gmx_atom_p[_id]
       _new_txt.append('%d  %s  %d  %s  %s  %d  %f  %f'%\
       (aid_b+_cnt,_atyp,\
        resid, resnm,\
        _anm, cgrp,\
        _par[1], _par[0])\
        +' ; added atom and borrowed from %s\n'%_atypB)
       print('Note: nonbonded parameters for atom type',_atyp,\
          'borrowed from type', _atypB)
       self.note+=1
      else:
     #  rl = self.atomtxt[_i]
       _new_txt.append('%d  %s  %d  %s  %s  %d  %f  %f'%\
       (aid_b+_cnt,_atyp,\
        resid, resnm,\
        _anm, cgrp,\
        0.0, guess_mass[_atyp[0]])\
        +' ; added atom but missing param !\n')
       print('Warning: missing parameters for atom type',_atyp)
       self.warn+=1

# insertion end
    if _i>=_natom: continue 
    rl = self.atomtxt[_i]
    _new_txt.append('%d'%\
       (self.aid[_i]+ _nadd)\
      +rl[find_end_idx_sec(rl, 1):])


    continue
 
   if _i>=_natom: continue
 
   _anm = self.atomnm[_i]
   _atyp, _atypB = mod_check.lookup(_anm)

   if _atyp is None:
    _new_txt.append(self.atomtxt[_i])
   elif _atyp== self.atyp[_i]:
    _new_txt.append(self.atomtxt[_i])
   elif _atyp in gmx_atom:
    rl = self.atomtxt[_i]
    _id = gmx_atom.index(_atyp)
    _par = gmx_atom_p[_id]
    _new_txt.append('%d  %s  %d  %s  %s  %d  %f  %f'%\
       (self.aid[_i],_atyp,\
        self.resid[_i], self.resnm[_i],\
        self.atomnm[_i], self.cgrp[_i],\
        _par[1], _par[0])\
      +rl[find_end_idx_sec(rl, 8):])
   elif _atypB in gmx_atom:
    rl = self.atomtxt[_i]
    _id = gmx_atom.index(_atypB)
    _par = gmx_atom_p[_id]
    _new_txt.append('%d  %s  %d  %s  %s  %d  %f  %f'%\
       (self.aid[_i],_atyp,\
        self.resid[_i], self.resnm[_i],\
        self.atomnm[_i], self.cgrp[_i],\
        _par[1], _par[0])\
      +' ; borrowed from %s\n'%_atypB)
    print('Note: nonbonded parameters for atom type',_atyp,\
          'borrowed from type', _atypB)
    self.note+=1
   else:
    rl = self.atomtxt[_i]
    _new_txt.append('%d  %s  %d  %s  %s  %d  %f  %f'%\
       (self.aid[_i],_atyp,\
        self.resid[_i], self.resnm[_i],\
        self.atomnm[_i], self.cgrp[_i],\
        0.0, guess_mass[_atyp[0]])\
      +' ; missing param !\n')
    print('Warning: missing parameters for atom type',_atyp)
    self.warn+=1
  self.atoms_update(  _new_txt   )

 def update_old_index(self,id_map,\
              my_ff, new_ff, atom_typ_changed=None):

  for _ih in range(len(self.sechead)):
   _sechead = self.sechead[_ih]
   _sec_st = self.secpos[_ih]
   _sec_end = self.secpos[_ih+1]
   if not _sechead in substitute_head: continue
   _nitem = substitute_nitem[substitute_head.index(_sechead)]

   for _ri in range(_sec_st, _sec_end, 1):
    _rl = self.txt[_ri]
    if _is_nondata_line(_rl): continue
    
# check only those non-commetn part of line
    _srl_nc = (_rl[:-1].split(';')[0]).split()

    if _sechead=='bonds':
     gmx_entry = my_ff.gmx_bonds
     gmx_param = my_ff.gmx_bond_params
    elif _sechead =='angles':
     gmx_entry = my_ff.gmx_angles
     gmx_param = my_ff.gmx_angle_params
    elif _sechead == 'pairs':
     gmx_entry = my_ff.gmx_pairs
     gmx_param = my_ff.gmx_pair_params
    elif _sechead == 'dihedrals':
     gmx_entry = my_ff.gmx_dihes
     gmx_param = my_ff.gmx_dihe_params
    elif _sechead == 'cmap':
     pass
    elif _sechead=='exclusions':
     pass
    else:
     continue


    _idx_entry = [int(_srl_nc[__i]) for __i in range(_nitem)]
    _rl_left = _rl[find_end_idx_sec(_rl, _nitem):]
    _mod_rl=''
    _atyp_entry = []
    _new_idx_entry = []
    atyp_nocheck = True
    for __j in range(_nitem):
     _mod_rl+='%d '%id_map[_idx_entry[__j]-1]
     _atyp_entry.append(self.atyp[id_map[_idx_entry[__j]-1]-1])
     _new_idx_entry.append(id_map[_idx_entry[__j]-1])

     if not atom_typ_changed is None:
      if id_map[_idx_entry[__j]-1] in atom_typ_changed:
       atyp_nocheck = False

    _mod_rl+=_rl_left
#    print(_mod_rl, end='')

    if _sechead=='exclusions'\
        or atyp_nocheck:
     self.txt[_ri] = _mod_rl
     continue

    if _sechead=='cmap':
     if atyp_nocheck:
      self.txt[_ri] = _mod_rl
      continue
     else:
      self.txt[_ri] = _mod_rl[:-1] + ' ; missing cmap param for '+\
         str[_atyp_entry]+'\n'
      print('Warning: cmap param for type',_atyp_entry,\
           _atyp_entry,' with missing parameters.')
      self.warn+=1
      continue

    if _sechead=='dihedrals':
     p_list = readFF.match_FF_list([_atyp_entry[1],_atyp_entry[2]], gmx_entry)
    else:
     p_list = readFF.match_FF_list(_atyp_entry, gmx_entry)

    if len(p_list)>0:
     if len(_srl_nc)> _nitem+1:
      _mod_rl= _mod_rl[:-1]+\
        ' ; new types '+str(_atyp_entry)+' defined but overriden here\n'
      print('Note: term',_sechead,'for type',_atyp_entry,_new_idx_entry,\
      'defined in FF but overriden with old params specified in topology')
      self.note+=1
     else:
      _mod_rl= _mod_rl[:-1]+\
        ' ; new types '+str(_atyp_entry)+'\n'

    elif len(_srl_nc)> _nitem+1:
      _mod_rl= _mod_rl[:-1]+\
        ' ; new types '+str(_atyp_entry)+' NOT defined but defined here\n'
      print('Note: term',_sechead,'for type',_atyp_entry,_new_idx_entry,\
      'not defined in FF but in topology with old parameters')
      self.note+=1
    else:
      _mod_rl= _mod_rl[:-1]+\
        ' ; new types '+str(_atyp_entry)+' NOT defined!\n'
      print('Warning: term',_sechead,'for type',_atyp_entry,_new_idx_entry,\
      'NOT defined in FF or in topology')
      self.warn+=1
    self.txt[_ri] = _mod_rl

 def insert_new_improper(self,conn,_add_atom_idx_list,xp):
  _is_add = np.zeros(len(self.atomnm), dtype=np.int32)
  _is_add[_add_atom_idx_list] = 1
  _secidx = self.sechead.index('dihedrals')
  _insert_txt = []
  atyp = self.atyp
  _insert_txt.append(';improper\n')
  for _i in range(len(self.atomnm)):

   _conn = conn[_i]
   if len(_conn)==3:
     if _is_add[[_i,_conn[0],_conn[1],_conn[2]]].sum()==0: continue
     if _is_add[_i]==0:
      print('Note: atom',_i+1, self.atomnm[_i],\
        'has three neighbors some of which are newly added')
      print('Note: this means its angle terms and imporper terms may not be well defined') 
     _dih = geom.dihedral(xp[_i],xp[_conn[0]],xp[_conn[1]],xp[_conn[2]])
     _k = 150.0 if not readTop.word_in_line(atyp[_i], aro_nm) is None else 300.0
     _insert_txt.append('%d   %d  %d  %d  2   %f  %f ; added imporper \n'%\
     (_i+1,_conn[0]+1,_conn[1]+1,_conn[2]+1,\
     _dih,_k))


  self.insert_txt_to_section(_insert_txt, _secidx)

 def insert_new_top(self,conn,_add_atom_idx_list, my_ff,xp,\
          dihStr=None, pairStr=None, debug=False):
# parameters for bonds : r0 from input geometry, k0 125000 default for PACE
# angles: theta0 from input, k0 300.0 kJ.mol defalt for PACE
# pairs, dihedrals: need to read from FF or added by users<--this part need to be paid attention to by users
# dihedrals type 2: 300.0 normally, 150.0 for aromatic
# dih parameters loaded from dihStr instead of connectivity
# and FF with dihStr is not None. This is a preferred way to generate
# topology 

  #print(dihStr)
   
  _is_add = np.zeros(len(self.atomnm), dtype=np.int32)
  _is_add[_add_atom_idx_list] = 1

  _secidx = self.sechead.index('bonds')
  _insert_txt = []
  for _i in conn.c12:
    if _is_add[_i].sum()>0: 
     _insert_txt.append('%d    %d  1   %f  %f ; added bond \n'%\
      (_i[0]+1,_i[1]+1,geom.bond(xp[_i[0]],xp[_i[1]])*0.1,125000.0))

  self.insert_txt_to_section(_insert_txt, _secidx)


  



  _secidx = self.sechead.index('pairs')
  _insert_txt = []

  if not pairStr is None:
     _insert_txt = [__j for __j in pairStr]
#     self.insert_txt_to_section(_insert_txt, _secidx)
     if 'exclusions' in self.sechead:
       _secidx_exc = self.sechead.index('exclusions')
       _insert_exc_txt = []
       for __k in  pairStr:
         srl = __k[:-1].split()
         _insert_exc_txt.append('%d  %d  ; added pair\n'\
              %(int(srl[0]),int(srl[1])))
       self.insert_txt_to_section(_insert_exc_txt, _secidx_exc)
     else:
       self.txt.append('[ exclusions ]\n')
       for __k in  pairStr:
            srl = __k[:-1].split()
            self.txt.append('%d  %d  ; added pair\n'\
              %(int(srl[0]),int(srl[1])))
       self.sechead.append('exclusions')
       _l_pos = self.secpos[-1]
       self.secpos.append(_l_pos+len(pairStr)+1) 
#  else:

  atyp = self.atyp
  _ff_pair = my_ff.gmx_pairs
  _ff_pair_param = my_ff.gmx_pair_params
  _store_pair=[]
  for _i in conn.c14:
     if _is_add[_i].sum()==0: continue
     _atyp_pair = [atyp[_i[0]], atyp[_i[3]]]
     _atyp_pair_r = [atyp[_i[3]], atyp[_i[0]]]

     if [_i[0],_i[3]] in _store_pair or\
        [_i[3],_i[0]] in _store_pair:
       continue
     else:
      _store_pair.append([_i[0],_i[3]])
      _store_pair.append([_i[3],_i[0]])

     _pair_id = -1
     if _atyp_pair in _ff_pair:
       _pair_id = _ff_pair.index(_atyp_pair)
     elif _atyp_pair_r in _ff_pair:
       _pair_id = _ff_pair.index(_atyp_pair_r)

     if _pair_id<0:    
      _insert_txt.append('%d   %d  1 ; added pair %s  %s missing param\n'%\
      (_i[0]+1, _i[3]+1, atyp[_i[0]], atyp[_i[3]]))
      print('Warning: added pair',_atyp_pair,_i[0]+1,_i[3]+1,\
           'with missing parameters')
      self.warn+=1
     else:
      _insert_txt.append('%d   %d  1   %g  %g ; added pair %s  %s\n'%\
      (_i[0]+1, _i[3]+1,\
       _ff_pair_param[_pair_id][0],_ff_pair_param[_pair_id][1],\
        atyp[_i[0]], atyp[_i[3]]))

  self.insert_txt_to_section(_insert_txt, _secidx)

  _secidx = self.sechead.index('angles')
  _insert_txt = []
  for _i in conn.c13:
    if _is_add[_i].sum()>0:
      _insert_txt.append('%d    %d  %d  1   %f  %f ; added angle\n'%\
      (_i[0]+1,_i[1]+1,_i[2]+1,\
             geom.angle(xp[_i[0]],xp[_i[1]],xp[_i[2]]),300.0))
  self.insert_txt_to_section(_insert_txt, _secidx)

  _secidx = self.sechead.index('dihedrals')
  _insert_txt = []
  _ff_dihe = my_ff.gmx_dihes
  _ff_dihe_param = my_ff.gmx_dihe_params

  if not dihStr is None:
   _insert_txt = [__j for __j in dihStr]
   #print(_insert_txt)
   self.insert_txt_to_section(_insert_txt, _secidx)
   return


  for _i in conn.c14:
   if _is_add[_i].sum()==0: continue
   _atyp_dihe = [atyp[_i[1]], atyp[_i[2]]]
   _atyp_dihe_r = [atyp[_i[2]], atyp[_i[1]]]

   _dihe_id = -1
   if _atyp_dihe in _ff_dihe:
     _dihe_id = _ff_dihe.index(_atyp_dihe)
   elif _atyp_dihe_r in _ff_dihe:
     _dihe_id = _ff_dihe.index(_atyp_dihe_r)

   if _dihe_id<0:
    _insert_txt.append('%d   %d  %d  %d  1 ;added dihedral %s  %s with missing param\n'%\
    (_i[0]+1,_i[1]+1,_i[2]+1, _i[3]+1, atyp[_i[1]], atyp[_i[2]]))
    print('Warning: added dihedral type',_atyp_dihe,\
       _i[0]+1,_i[1]+1,_i[2]+1, _i[3]+1,'with missing parameters')
    self.warn+=1
   else:
    _insert_txt.append('%d   %d  %d  %d  1   %f  %f  %d  ; added dihedral %s  %s\n'%\
    (_i[0]+1, _i[1]+1, _i[2]+1, _i[3]+1,\
         _ff_dihe_param[_dihe_id][0],_ff_dihe_param[_dihe_id][1],\
      _ff_dihe_param[_dihe_id][2],\
      atyp[_i[1]], atyp[_i[2]]))


  self.insert_txt_to_section(_insert_txt, _secidx)



 def insert_txt_to_section(self, _insert_txt, _secidx):
  _new_txt = self.txt[:self.secpos[_secidx+1]]
  for _rl in _insert_txt:
   _new_txt.append(_rl)
 
  _new_txt+= self.txt[self.secpos[_secidx+1]:]
  _nadd = len(_insert_txt)
  for _is in range(_secidx+1, len(self.secpos),1):
   self.secpos[_is]+= _nadd
  
  self.txt = _new_txt

 def show_section(self,sechead):
  _secidx = self.sechead.index(sechead)
  for _rl in self.txt[self.secpos[_secidx]:self.secpos[_secidx+1]]:
   print(_rl,end='')
 
 def output(self, fnm=None):

  
  _atom_head_idx = self.sechead.index('atoms')
#  _comp_txt = self.txt[:self.secpos[_atom_head_idx]]+['[ atoms ]\n']\
#     +self.atomtxt+\
#      self.txt[self.secpos[_atom_head_idx+1]:]  

  _comp_txt = self.txt[:self.secpos[_atom_head_idx]]+['[ atoms ]\n']\
     +self.atomtxt

  _sec_order = ['bonds','exclusions','pairs','angles',\
            'dihedrals','cmap','position_restraints','system','molecules']

  for _is in _sec_order:
   if _is in self.sechead:
    for __k in range(len(self.sechead)):
     if self.sechead[__k]==_is:
      
      _comp_txt+= self.txt[self.secpos[__k]:self.secpos[__k+1]]

  

  if fnm is None:
   return _comp_txt
  else:
   f=open(fnm,'w')
   for _rl in _comp_txt:
    f.write(_rl)
   f.close() 
