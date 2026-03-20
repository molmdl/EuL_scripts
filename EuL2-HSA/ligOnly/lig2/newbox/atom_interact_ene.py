import sys
import numpy as np

def load_param(ftops,fparam,pace_ff=False):

# ftops inlcudes names of a list of top files, each
# containing atomtyping info for each mol typ
# the order of these files in the list should
# be consistent with the order of mol typs in trj
# fparam includes a list of names of param files

 typnm=[]
 _mol_size=[]
 charge=[]
 _n_top = len(ftops)

 for _i in range(_n_top):
  _fnm = ftops[_i]

  _tn=[]
  _rd = False
  _ct=[]
  for rl in open(_fnm,'r'):
   if 'atoms' in rl and '[' in rl:
    _rd = True
    continue

   srl=rl[:-1].split()

   if len(srl)==0: continue
   if ';' in srl[0]: continue

   if '[' in srl[0] and _rd: break

   if _rd:
    _tn.append(srl[1])
    _ct.append(float(srl[6]))

  charge.append(_ct)
  typnm.append(_tn)
  _mol_size.append(len(_tn))

#look up corresponding vdw param from fparam file

 _nm = []
 _v = []
 _w = []


 for _f in  fparam:
  _rd = False
  for rl in open(_f,'r'):
   if '[' in rl and 'atomtype' in rl:
    _rd = True
    continue

   srl=rl[:-1].split()

   if len(srl)==0: continue
   if ';' in srl[0]: continue
# macro in ff file will be ignored
# but this could lead to wrong
# interpretation of FF. Caution needed!
   if '#' in srl[0]: continue
   if '[' in srl[0] and _rd: break

   if _rd:

    _nm.append(srl[0])
    if pace_ff:
     # C6
     _v.append(float(srl[4]))
     # c12
     _w.append(float(srl[5]))

    else:
     _v.append(float(srl[5]))
     _w.append(float(srl[6]))

# load nb matrix
 _n_tp = len(_nm)
 _is_mat = np.zeros((_n_tp, _n_tp), dtype=np.int16)
 _v_mat =  np.zeros((_n_tp, _n_tp))
 _w_mat =  np.zeros((_n_tp, _n_tp))

 for _f in  fparam:
  _rd = False
  for rl in open(_f,'r'):
   if '[' in rl and 'nonbond-param' in rl:
    _rd = True
    continue

   srl=rl[:-1].split()

   if len(srl)==0: continue
   if ';' in srl[0]: continue
# macro in ff file will be ignored
# but this could lead to wrong
# interpretation of FF. Caution needed!
   if '#' in srl[0]: continue
   if '[' in srl[0] and _rd: break

   if _rd:

    _idxa = _nm.index(srl[0])
    _idxb = _nm.index(srl[1])
    _is_mat[_idxa,_idxb] =1
    _is_mat[_idxb,_idxa] =1
    _v_mat[_idxa,_idxb] = float(srl[3])
    _w_mat[_idxa,_idxb] = float(srl[4])
    _v_mat[_idxb,_idxa] = _v_mat[_idxa,_idxb]
    _w_mat[_idxb,_idxa] = _w_mat[_idxa,_idxb]
    


 _vdw_v=[]
 _vdw_w=[]
 _nb_mat = []
 _nb_v_mat = []
 _nb_w_mat = []
# print(_n_tp)
 for _i in range(_n_top):
#  print(_i)
  _tn = typnm[_i]

  _vt = []
  _wt = []


  for _j in _tn:
   if not _j in _nm:
    print( 'Error. atomtyp', _j,'in mol typ',_i,'not found in FF param database')
    sys.exit()

   _idx = _nm.index(_j)
   _vt.append(_v[_idx])
   _wt.append(_w[_idx])

  _vdw_v.append(_vt)
  _vdw_w.append(_wt)    


 for _i in range(_n_top):

  _tn = typnm[_i]


  _nbt = []
  _nbvt=[]
  _nbwt=[]


  for _j in range(_n_top):
   _tm = typnm[_j]
   _ni = len(_tn)
   _nj = len(_tm)
   _c_nb = np.zeros((_ni, _nj),dtype=np.int16)
   _c_v = np.zeros((_ni, _nj))
   _c_w = np.zeros((_ni, _nj))

   for _ki,_k in enumerate( _tn):
    for _lj,_l in enumerate( _tm):
     _idxi = _nm.index(_k)
     _idxj = _nm.index(_l)
     _c_nb[_ki,_lj] = _is_mat[_idxi, _idxj]
     _c_v[_ki,_lj] = _v_mat[_idxi, _idxj]
     _c_w[_ki,_lj] = _w_mat[_idxi, _idxj]

   _nbt.append(_c_nb)
   _nbvt.append(_c_v)
   _nbwt.append(_c_w)

  _nb_mat.append(_nbt)
  _nb_v_mat.append(_nbvt)
  _nb_w_mat.append(_nbwt)

 return (_mol_size, charge, _vdw_v, _vdw_w, _nb_mat, _nb_v_mat, _nb_w_mat)


def parameter_mat(_mol_size, charge, vdw_v, vdw_w, nb_mat, nb_v_mat, nb_w_mat,\
                   _ele_f = 138.935, LJrule =2, ):
# load parameters taken from FF and convert it into param matrix
# for col and vdw and vdw in C6 and C12 format
# how to convert params depends on LJrule
# 1 for pace ; 2 for charmm
 _n_tp=len(_mol_size)
 if _n_tp!=len(charge)\
 or _n_tp!=len(vdw_v)\
 or _n_tp!=len(vdw_w):
  print( 'Error! unmatched number of types')
  sys.exit()

 for _i in range(_n_tp):
  if _mol_size[_i]!=len(charge[_i])\
 or  _mol_size[_i]!=len(vdw_v[_i])\
 or _mol_size[_i]!=len(vdw_w[_i]):
   print( 'Error! unmatched atom count for typ',_i)
   sys.exit()

 col_param_mat=[]
 vdw_param_C6_mat=[]
 vdw_param_C12_mat=[]


 for _i in range(_n_tp):
  _c=[]
  _v6=[]
  _v12=[]

  for _j in range(_n_tp):

   _ct= np.zeros((_mol_size[_i],_mol_size[_j]),dtype=np.float64)
   _v6t= np.zeros((_mol_size[_i],_mol_size[_j]),dtype=np.float64)
   _v12t= np.zeros((_mol_size[_i],_mol_size[_j]),dtype=np.float64)

   _nb = nb_mat[_i][_j]
   _nb_v = nb_v_mat[_i][_j]
   _nb_w = nb_w_mat[_i][_j]     

   _ch_i = charge[_i]
   _ch_j = charge[_j]
   _vv_i = vdw_v[_i]
   _vv_j = vdw_v[_j]
   _vw_i = vdw_w[_i]
   _vw_j = vdw_w[_j]



   for _ik in range(_mol_size[_i]):

    for _jk in range(_mol_size[_j]):

     _ct[_ik,_jk] = _ele_f * _ch_i[_ik]*_ch_j[_jk]



     if LJrule==2:
      if _nb[_ik, _jk]==1:
       _v_v = _nb_v[_ik, _jk]
       _v_w = _nb_w[_ik, _jk]
      else:
       _v_v = (_vv_i[_ik]+_vv_j[_jk])*0.5
       _v_w = (_vw_i[_ik]*_vw_j[_jk])**0.5

      _v6t[_ik,_jk] = 4.*_v_w*_v_v**6.0
      _v12t[_ik, _jk] = 4.*_v_w*_v_v**12.0
     elif LJrule==1:
      if _nb[_ik, _jk]==1:
       _v6t[_ik,_jk] = _nb_v[_ik, _jk]
       _v12t[_ik,_jk] = _nb_w[_ik, _jk]
      else:
       _v6t[_ik,_jk] = (_vv_i[_ik]*_vv_j[_jk])**0.5
       _v12t[_ik,_jk] = (_vw_i[_ik]*_vw_j[_jk])**0.5
     else:
# nb mat dont work for other LJrules
      _C6_i = 4.*_vw_i[_ik]*_vv_i[_ik]**6.0
      _C6_j = 4.*_vw_j[_jk]*_vv_j[_jk]**6.0
      _C12_i = 4.*_vw_i[_ik]*_vv_i[_ik]**12.0
      _C12_j = 4.*_vw_j[_jk]*_vv_j[_jk]**12.0

      _v6t[_ik,_jk] = (_C6_i*_C6_j)**0.5
      _v12t[_ik, _jk] = (_C12_i*_C12_j)**0.5

   _c.append(_ct)
   _v6.append(_v6t)
   _v12.append(_v12t)

     
  col_param_mat.append(_c)
  vdw_param_C6_mat.append(_v6)
  vdw_param_C12_mat.append(_v12)



 return (col_param_mat, vdw_param_C6_mat, vdw_param_C12_mat) 

def gen_RF_param (_cut_off, e_rf):
 return [(e_rf-1.0)/(2*e_rf+1.0)/(_cut_off**3.0),\
         3*e_rf/(2*e_rf+1.0)/_cut_off]


def atom_inter_ene (_mol_cnt,_mol_size, _ri_1,_ri_6,_ri_12,\
                    _param_1, _param_6, _param_12, reactionField=False,\
                    RF_param=None, _r2 = None, _contact=None, byMol=False):

 if reactionField:
  RF_p1= RF_param[0]
  RF_p2= RF_param[1]

 _n_tp=len(_mol_size)
 if byMol:
  _n_mol=sum(_mol_cnt)
  _m_st, _m_tp=set_mol_atom_idx(_mol_cnt,_mol_size) 

 col_inter=[]
 vdw_inter=[]

 for _ii in range(_n_mol if byMol else _n_tp):
  _i = _m_tp[_ii] if byMol else _ii
  _c=[]
  _v=[]
  for _j in range(_n_tp):
   if reactionField:
    _c.append((_ri_1[_ii][_j] + RF_p1*_r2[_i][_j] - RF_p2*_contact[_ii][_j])*_param_1[_i][_j])
   else:
    _c.append(_ri_1[_ii][_j]*_param_1[_i][_j])
   _v.append(_ri_12[_ii][_j]*_param_12[_i][_j] - _ri_6[_ii][_j]*_param_6[_i][_j])

  col_inter.append(_c)
  vdw_inter.append(_v)

 return (col_inter, vdw_inter)
 


def frag_frag_distance (_xp1, xp2, _mark, edge_len=[1e+6,1e+6,1e+6], _cut = 4.5,\
                         onlyContact=True, singleSite = False, dualBoundary=False,\
                         _cut_sec= 4.5, returnDx=False):
# cal contact matrix for coord of fragment 1 (xp1) and fragment 2 (xp2)
# _mark: flag to mark which atom of frag 2 will be investigated when checking if it is in vicinity of frag 1
# onlyContact: ture only contact matrix false: all distance matrix needed for energy evalulation  

# available option by flags
# onlyContact = True and singleSite =True: input xp1 is a verctori for reference point xp1 not a list of vectors
#                   return adjacency vector, distance square vector, a list of dx for other vectors in xp2 
# onlyContact = False and singleSite = True: return adjecency vectorm distance square vector, distance vector,
#                                             r^6 vector, r^12 vector
# onlyContact = True and singeSite = False and dualBoundary = False: xp1 is a list of vector
#                                   if returnDx=True: return adjacency matrix, distance square matrix, matrix of dx vectors
#                                    else: adjacency matrix
# onlyContact = False and singeSite = False: xp1 is a list of vector
#                                          return adjacency matrix, distance suqare matix, distance matrix, r^6 matrix, r^12 matrix
# onlyContact = True and singeSite = False and dualBoundary = True: xp1 is a list of vector
#                                   if returnDx=true: return adjacency matrix for cut1, adjacency matrix for cut2
#                                                            distance square matrix for cut1, matrix of dx vectors for cut1
#                                           else: adjacency matrix for cut1, adjacency matrix for cut2 

 if singleSite:
  xp1 = np.zeros((1,3))
  xp1[0,:] = _xp1[:]
  _n_x1 = 1
 else:
  xp1=_xp1
  _n_x1 = len(xp1)

 _n_x2 = len(xp2)
 _cut2 = _cut*_cut
 _cut2_sec = _cut_sec*_cut_sec 

 _r2 = np.zeros(_n_x2, dtype=np.float64)

 if not onlyContact: 
  _con_stat_1 = np.zeros((_n_x1, _n_x2),dtype=np.float64)
  _con_stat_6 = np.zeros((_n_x1, _n_x2),dtype=np.float64)
  _con_stat_12 = np.zeros((_n_x1, _n_x2),dtype=np.float64)


 _con_stat2 = np.zeros((_n_x1, _n_x2),dtype=np.float64)

 _con_stat = np.zeros((_n_x1, _n_x2),dtype=np.float64)

 if dualBoundary:
  _con_stat_sec = np.zeros((_n_x1, _n_x2),dtype=np.float64)


 if returnDx:
  _dx_all = np.zeros((_n_x1,_n_x2,3))

 for _j in range(_n_x1):

    if returnDx:
     _dx_from_i  = _dx_all[_j]
     _dx_from_i[:,:]  = xp2 -xp1[_j,:]
    else:
     _dx_from_i  =  xp2 -xp1[_j,:]

    _r2[:]=0.0

    for _d in range(3):
      _dx_t = _dx_from_i[:,_d]

      _idx_gt_hb = np.where(_dx_t>\
                                0.5*edge_len[_d])
      _dx_t[_idx_gt_hb]-= edge_len[_d]
      _idx_gt_hb = np.where(_dx_t<-0.5*\
                                edge_len[_d])
      _dx_t[_idx_gt_hb]+= edge_len[_d]

      _dx_from_i[:,_d] = _dx_t[:]

      _r2+= (_dx_t*_dx_t)

    if not onlyContact:
     _cs_1 = _con_stat_1[_j,:]
     _cs_6 = _con_stat_6[_j,:]
     _cs_12 = _con_stat_12[_j,:]



    _cs2 = _con_stat2[_j,:]
    _cs = _con_stat[_j,:]

    if dualBoundary:
     _cs_sec = _con_stat_sec[_j,:]

#    _idx_within_rc = np.where(np.absolute(_r2-(_cut2-0.2))<=0.2)

    _idx_within_rc = np.where(_r2<=_cut2)
    
    if dualBoundary:
     _idx_within_rc_sec = np.where(_r2<=_cut2_sec)

    _r2_cut = _r2[_idx_within_rc]+ 1e-8
    if not onlyContact:
     _r2_i = np.reciprocal(_r2_cut)
     _r1_i = np.sqrt(_r2_i)
#    _r1 = np.reciprocal(_r1_i)
     _r6_i = _r2_i*_r2_i*_r2_i
     _r12_i = _r6_i*_r6_i

     _cs_1[_idx_within_rc]+=(_r1_i*_mark[_idx_within_rc])
     _cs_6[_idx_within_rc]+=(_r6_i*_mark[_idx_within_rc])
     _cs_12[_idx_within_rc]+=(_r12_i*_mark[_idx_within_rc])



    _cs[_idx_within_rc]+=(_mark[_idx_within_rc])
    _cs2[_idx_within_rc]+=(_r2_cut*_mark[_idx_within_rc])

    if dualBoundary:
     _cs_sec[_idx_within_rc_sec]+=(_mark[_idx_within_rc_sec])

 if onlyContact:
  if singleSite:
   return _con_stat[0],_con_stat2[0],_dx_from_i
  else:
   if dualBoundary:
    if returnDx:
     return (_con_stat, _con_stat_sec, _con_stat2, _dx_all)
    else:
     return (_con_stat, _con_stat_sec)
   else:
    if returnDx:
     return (_con_stat,_con_stat2, _dx_all)
    else:
     return _con_stat
 else:
  if singleSite:
   return (_con_stat[0],_con_stat2[0], _con_stat_1[0], _con_stat_6[0],_con_stat_12[0])
  else:
   return (_con_stat,_con_stat2, _con_stat_1, _con_stat_6,_con_stat_12)


def set_mol_atom_idx(_mol_cnt,_mol_size):
 _n_tp = len(_mol_cnt)
# print _mol_cnt_cal
 __st_m =0
 _m_st=[]
 _m_tp=[]
 for _i in range(_n_tp):
#starting atom idx for each typ of mol

  for _j in range(_mol_cnt[_i]):

   _m_tp.append(_i)
   _m_st.append(__st_m)
   __st_m+=_mol_size[_i]

#starting atom idx for each mol  
 _m_st.append(__st_m)
 return (_m_st, _m_tp)

def init_atom_mat_by_mol_typ(_mol_cnt,_mol_size):
 _re = []
 _n_tp = len(_mol_cnt)

 
 for _i in range(_n_tp):
  _s = _mol_size[_i]
  for _j in range(_mol_cnt[_i]):
   
   _re.append([np.zeros((_s, _k)) for _k in _mol_size ])

 return _re

def dist_per_mol_for_typ(xp, _mol_cnt,_mol_size, edge_len=[1e+6,1e+6,1e+6], _cut = 4.5,\
                           reactionField=False):

 _con_stat,_con_stat2, _con_stat_1, _con_stat_6,_con_stat_12 =\
      frag_frag_distance (xp, xp, np.ones(len(xp),dtype=np.int16),\
                         edge_len=edge_len, _cut = _cut,\
                         onlyContact=False)  
 _m_st, _m_tp = set_mol_atom_idx(_mol_cnt,_mol_size)
 _n_mol = sum(_mol_cnt)
 _n_tp = len(_mol_size)
# for _i,_j in enumerate(_con_stat): print(_i, _j.sum())
# print(_m_st, _m_tp)
 
# output will be a series of atomic contact matrices betwen a mol with its surrounds
# for mol i (with a atoms) and there are N types of mols and Xth type of mols have b atom
# then matrix re[i][X][a,b] denotes atomic interaction matrix between mol i and all atoms of type X
 re = init_atom_mat_by_mol_typ(_mol_cnt,_mol_size)
 re2 = init_atom_mat_by_mol_typ(_mol_cnt,_mol_size)
 re_1 = init_atom_mat_by_mol_typ(_mol_cnt,_mol_size)
 re_6 = init_atom_mat_by_mol_typ(_mol_cnt,_mol_size)
 re_12 = init_atom_mat_by_mol_typ(_mol_cnt,_mol_size)

 for _i in range(_n_mol):
  _i_st = _m_st[_i]
  _i_end=  _m_st[_i+1]
  for _j in range(_n_mol):
   if _i==_j: continue
   _j_tp = _m_tp[_j]
   _j_st = _m_st[_j]
   _j_end=  _m_st[_j+1]

   re[_i][ _j_tp]+= _con_stat[_i_st: _i_end, _j_st: _j_end]
   re2[_i][ _j_tp]+= _con_stat2[_i_st: _i_end, _j_st: _j_end]
   re_1[_i][ _j_tp]+= _con_stat_1[_i_st: _i_end, _j_st: _j_end]
   re_6[_i][ _j_tp]+= _con_stat_6[_i_st: _i_end, _j_st: _j_end]
   re_12[_i][ _j_tp]+= _con_stat_12[_i_st: _i_end, _j_st: _j_end]

 if reactionField:
  return (re,re2, re_1, re_6,re_12)
 else:
  return (re, re_1, re_6,re_12)



def atom_distance(xp, _mol_cnt,_mol_size,_to_cal,edge_len=[1e+6,1e+6,1e+6], _cut = 4.5,\
                  reactionField=False, _excludeList=None):

# Calculate average 1/r, 1/r^6, and 1/r^12 between atom pairs of
# mol typ A and mol typ B
#
# out: martices
# M_AB(i,j) for each mol typ A
# M_AB(i,j)= SIGMA_J=1..NB{1/R_A_i-B_J_j} 
# R_A_i-B_J_j funtion between i atom in a given mol of typ A and j atom in Jth mol of typ B 
# R can be r, r^6, r^12
# M_AB(i,j) average over all given mols of typ A

 _n_x = len(xp)
 _cut2 = _cut*_cut

 _r2 = np.zeros(_n_x, dtype=np.float64)
 _mark = np.zeros(_n_x, dtype=np.int32)
#used to calculate col, vdw 6 and 12 terms
 _re_1=[]
 _re_6=[]
 _re_12=[]
 _re=[]
 _re2=[]

 _n_tp = len(_mol_cnt)

 _mol_st = np.zeros(_n_tp, dtype=np.int32)
#num of mols of a typ that will be considered for calculation 
 _mol_cnt_cal = np.zeros(_n_tp, dtype=np.int32)

 _n_mol = sum(_mol_cnt)
 _mol_tp_id = np.zeros(_n_mol, dtype=np.int32) 
 _mark_mol= np.zeros(_n_mol, dtype=np.int32)
 _s = 0
 for _i in range(_n_tp):
#typ id for each mol: mols of same type are continuous in order
  _mol_tp_id[_s:_s+_mol_cnt[_i]]=_i
  _s+= _mol_cnt[_i]

 for _i in _to_cal: _mol_cnt_cal[_mol_tp_id[_i]]+=1
# print _mol_cnt_cal
 __st_m =0
 _m_st=[]

 __st=0
 for _i in range(_n_tp):
#starting atom idx for each typ of mol
  _mol_st[_i]=__st
  __st+=_mol_cnt[_i]*_mol_size[_i]

  for _j in range(_mol_cnt[_i]):


   _m_st.append(__st_m)
   __st_m+=_mol_size[_i]

#starting atom idx for each mol  
 _m_st.append(__st_m)
# print _m_st
# to mark each atom if it is considered for interaction calculation
 _mark_cal = np.zeros(_n_x, dtype=np.int32)
 for _i in _to_cal:
  _mark_cal[_m_st[_i]:_m_st[_i+1]]=1
  _mark_mol[_i]=1

 if not _excludeList is None:
  for _i in _excludeList:
   _mark_mol[_i]=0

 for _i in range(_n_tp):
  _rr=[]
  _rr2 =[]
  _rr_1=[]
  _rr_6=[]
  _rr_12=[]
  for _j in range(_n_tp):
   _rr.append(np.zeros((_mol_size[_i], _mol_size[_j]), dtype=np.float64))
   _rr2.append(np.zeros((_mol_size[_i], _mol_size[_j]), dtype=np.float64))
   _rr_1.append(np.zeros((_mol_size[_i], _mol_size[_j]), dtype=np.float64))
   _rr_6.append(np.zeros((_mol_size[_i], _mol_size[_j]), dtype=np.float64))
   _rr_12.append(np.zeros((_mol_size[_i], _mol_size[_j]), dtype=np.float64))

#   print _rr[-1].shape
  _re.append(_rr)
  _re2.append(_rr2)
  _re_1.append(_rr_1)
  _re_6.append(_rr_6)
  _re_12.append(_rr_12)

# print _mol_st
 _mol_idx=0
 _idx = 0
 for _ti in range(_n_tp):
#distance matrix between a mol of typ _ti with all other atoms
  _con_stat_1 = np.zeros((_mol_size[_ti], _n_x),dtype=np.float64)
  _con_stat_6 = np.zeros((_mol_size[_ti], _n_x),dtype=np.float64)
  _con_stat_12 = np.zeros((_mol_size[_ti], _n_x),dtype=np.float64)
  _con_stat2 = np.zeros((_mol_size[_ti], _n_x),dtype=np.float64)
  _con_stat = np.zeros((_mol_size[_ti], _n_x),dtype=np.float64)
  for _i in range(_mol_cnt[_ti]):
# check if need to be cal
   if _mark_mol[_mol_idx]==0:
    _mol_idx+=1
    _idx+= _mol_size[_ti]
    continue

   _mol_idx+=1

   _mark[:]=1
#strating atom idx of this mol
   _st = _mol_st[_ti] + _mol_size[_ti]*_i
#used to exlude self-interactions
   _mark[_st:_st + _mol_size[_ti]]=0

   for _j in range(_mol_size[_ti]):


    _dx_from_i  =  xp -xp[_idx,:]

    _r2[:]=0.0

    for _d in range(3):
      _dx_t = _dx_from_i[:,_d]

      _idx_gt_hb = np.where(_dx_t>\
                                0.5*edge_len[_d])
      _dx_t[_idx_gt_hb]-= edge_len[_d]
      _idx_gt_hb = np.where(_dx_t<-0.5*\
                                edge_len[_d])
      _dx_t[_idx_gt_hb]+= edge_len[_d]


      _r2+= (_dx_t*_dx_t)

    _cs2 = _con_stat2[_j,:]
    _cs_1 = _con_stat_1[_j,:]
    _cs_6 = _con_stat_6[_j,:]
    _cs_12 = _con_stat_12[_j,:]
    _cs = _con_stat[_j,:]

#    _idx_within_rc = np.where(np.absolute(_r2-(_cut2-0.2))<=0.2)

    _idx_within_rc = np.where(_r2<=_cut2)
    _r2_cut = _r2[_idx_within_rc]+ 1e-8
    _r2_i = np.reciprocal(_r2_cut)
    _r1_i = np.sqrt(_r2_i)
#    _r1 = np.reciprocal(_r1_i)
    _r6_i = _r2_i*_r2_i*_r2_i
    _r12_i = _r6_i*_r6_i

    _cs2[_idx_within_rc]+=(_r2_cut*_mark[_idx_within_rc]*_mark_cal[_idx_within_rc])    
    _cs_1[_idx_within_rc]+=(_r1_i*_mark[_idx_within_rc]*_mark_cal[_idx_within_rc])
    _cs_6[_idx_within_rc]+=(_r6_i*_mark[_idx_within_rc]*_mark_cal[_idx_within_rc])
    _cs_12[_idx_within_rc]+=(_r12_i*_mark[_idx_within_rc]*_mark_cal[_idx_within_rc])
    _cs[_idx_within_rc]+=(_mark[_idx_within_rc]*_mark_cal[_idx_within_rc])

    _idx+=1

   
  for _tj in range(_n_tp):
   _r_c_1 = _re_1[_ti][_tj]
   _r_c_6 = _re_6[_ti][_tj]
   _r_c_12 = _re_12[_ti][_tj]
   _r_c = _re[_ti][_tj]
   _r_c2 = _re2[_ti][_tj]

#   print _r_c.shape
   _sta = _mol_st[_tj] 
   for _j in range(_mol_cnt[_tj]):
#    print _sta, _sta+_mol_size[_tj]
    _r_c_1+=_con_stat_1[:, _sta:_sta+ _mol_size[_tj]]
    _r_c_6+=_con_stat_6[:, _sta:_sta+ _mol_size[_tj]]
    _r_c_12+=_con_stat_12[:, _sta:_sta+ _mol_size[_tj]]
    _r_c+=_con_stat[:, _sta:_sta+ _mol_size[_tj]]
    _r_c2+=_con_stat2[:, _sta:_sta+ _mol_size[_tj]]

    _sta+= _mol_size[_tj]


 for _i in range(_n_tp):
  for _j in range(_n_tp):
   if _mol_cnt_cal[_i]==0 or _mol_cnt_cal[_j]==0:
    _re_1[_i][_j][:,:]=0.0
    _re_6[_i][_j][:,:]=0.0
    _re_12[_i][_j][:,:]=0.0
    _re[_i][_j][:,:]=0.0
    _re2[_i][_j][:,:]=0.0

    continue
   _re_1[_i][_j]*=(1./float(_mol_cnt_cal[_i]))
   _re_6[_i][_j]*=(1./float(_mol_cnt_cal[_i]))
   _re_12[_i][_j]*=(1./float(_mol_cnt_cal[_i]))
   _re[_i][_j]*=(1./float(_mol_cnt_cal[_i]))
   _re2[_i][_j]*=(1./float(_mol_cnt_cal[_i]))

 if not reactionField:
  return (_re, _re_1, _re_6,_re_12)
 else:
  return (_re,_re2, _re_1, _re_6,_re_12)



# test

if __name__=='__main__':


 xp = np.zeros((100,3))

 for _i in range(10):
  for _j in range(10):
   xp[_i*10+_j,0]= 5.0*_i
   xp[_i*10+_j,1]= 5.0*_j
   xp[_i*10+_j,2]= 0.0

 _mark = np.ones(100,dtype=np.int32)
 _mark[:3]=0

 _con, _con_sec = frag_frag_distance(xp[:3],xp,_mark,_cut=6.0,\
                edge_len=[50.0,50.0,50.0], dualBoundary=True, _cut_sec = 8.0)

 print(_con)

 for _c in _con:
  for _i in range(10):
   for _j in range(10):
    print(_c[_i*10+_j],end=' ')
   print()

  print()

 for _c in _con_sec:
  for _i in range(10):
   for _j in range(10):
    print(_c[_i*10+_j],end=' ')
   print()

  print()

