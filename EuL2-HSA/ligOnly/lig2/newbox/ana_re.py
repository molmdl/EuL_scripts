import numpy as np
import MDAnalysis
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import PDB
from atom_interact_ene import frag_frag_distance

def xp_box(xp, box):
# use the first atom as the ref
 _xp = xp + 0.0

 _dx = _xp - _xp[0]
 for d in range(3):
  _idx = np.where(_dx[:,d]>0.5*box[d])
  _xp[_idx, d]-= box[d]
  _idx = np.where(_dx[:,d]<-0.5*box[d])
  _xp[_idx, d]+= box[d]

 return _xp

def comass(xp):
 n_f = 1./float(len(xp))

 _com = np.zeros(3)


 for i in range(3):
  _com[i] = xp[:,i].sum()

 return _com*n_f


def find_idx( my_pdb,anm_list = None, start_=None, excl_ = None):
 _list=[]

 for _i,_j in enumerate(my_pdb.atomnms):
  _append = False
  if not anm_list is None:
   if _j in anm_list: _append=True

  if not start_ is None:
   if _j[0] == start_: _append = True

  if not excl_ is None:
   if _j==excl_: _append = False

  if _append: _list.append(_i)

 if _list==[]:
  print('Error: cannot find any atoms')
  sys.exit()


 return np.array(_list, dtype=np.int32)

def radial_dist_hist(xp_a,xp_b, box=[1e+6, 1e+6,1e+6], r_range = [0.5, 3.5] , dr=0.025):
# n_r = len(full_trj)
# n_m = len(full_trj[0])

 n_bin = int((r_range[1]-r_range[0])/dr)

 _rdf = np.zeros(n_bin)

# for _xp in full_trj:
 _adj, _dist2, _dist, _duma, _dumb =\
      frag_frag_distance(xp_a, xp_b, np.ones(len(xp_b), dtype=np.int16),\
          _cut = r_range[1], edge_len=box,\
          onlyContact =False, singleSite=False)
 np.fill_diagonal(_dist,0.0)
#   print(_xp)
#   print(_dist2)
#   sys.exit()
 _hist, _bine = np.histogram(np.sqrt(_dist2 + 1e-8), bins= n_bin, range=r_range)
 _rdf+=_hist

# print(_rdf)

 binc = 0.5*(_bine[1:]+_bine[:-1])
 _r_1 = np.reciprocal(binc)

 return(_rdf*_r_1*_r_1, binc)


##MAIN

xtc = 'v1.xtc' 
PDBname = 'v1.pdb' 

u=MDAnalysis.Universe(PDBname,xtc)
my_pdb = PDB.PDBdata(PDBname)

jump_step=1
Eu_id = find_idx(my_pdb,anm_list=['EU3'])
Np_id = find_idx(my_pdb,anm_list=['N6'])
N_id = find_idx(my_pdb, start_='N', excl_='N6')
O_id = find_idx(my_pdb, start_='O', excl_='OW')
OW_id = find_idx(my_pdb, anm_list = ['OW'])

first_=True
for ts in u.trajectory[::jump_step]:
    if ts.frame%100==0:
     print( 'time',u.trajectory.time,'fr id',ts.frame)
     sys.stdout.flush()

    xp = ts.positions
    edge_len = ts.dimensions[:3]

    if first_:
     _cnt=1
     _bin_Eu_N, _binc_Eu_N = radial_dist_hist(xp[Eu_id], xp[N_id], box=edge_len)
     _bin_Eu_O, _binc_Eu_O = radial_dist_hist(xp[Eu_id], xp[O_id], box=edge_len)
     _bin_Eu_OW, _binc_Eu_OW = radial_dist_hist(xp[Eu_id], xp[OW_id], box=edge_len)
     _bin_Eu_Np, _binc_Eu_Np = radial_dist_hist(xp[Eu_id], xp[Np_id], box=edge_len)
     first_=False
    else:

     _cnt+=1
     _bin_t , dummy = radial_dist_hist(xp[Eu_id], xp[N_id], box=edge_len)
     _bin_Eu_N+= _bin_t
     _bin_t , dummy = radial_dist_hist(xp[Eu_id], xp[O_id], box=edge_len)
     _bin_Eu_O+= _bin_t
     _bin_t , dummy = radial_dist_hist(xp[Eu_id], xp[OW_id], box=edge_len)
     _bin_Eu_OW+= _bin_t
     _bin_t , dummy = radial_dist_hist(xp[Eu_id], xp[Np_id], box=edge_len)
     _bin_Eu_Np+= _bin_t



print('Eu-N', (_bin_Eu_N*_binc_Eu_N).sum()/_bin_Eu_N.sum(), (_bin_Eu_N*_binc_Eu_N*_binc_Eu_N).sum()/_cnt)
print('Eu-O', (_bin_Eu_O*_binc_Eu_O).sum()/_bin_Eu_O.sum(), (_bin_Eu_O*_binc_Eu_O*_binc_Eu_O).sum()/_cnt)
print('Eu-OW', (_bin_Eu_OW*_binc_Eu_OW).sum()/_bin_Eu_OW.sum(), (_bin_Eu_OW*_binc_Eu_OW*_binc_Eu_OW).sum()/_cnt)
print('Eu-Npy', (_bin_Eu_Np*_binc_Eu_Np).sum()/_bin_Eu_Np.sum(), (_bin_Eu_Np*_binc_Eu_Np*_binc_Eu_Np).sum()/_cnt)

plt.figure()
plt.plot(_binc_Eu_N, _bin_Eu_N/_cnt,label='Eu-N')
plt.plot(_binc_Eu_O, _bin_Eu_O/_cnt,label='Eu-O')
plt.plot(_binc_Eu_OW, _bin_Eu_OW/_cnt,label='Eu-OW')
plt.plot(_binc_Eu_Np, _bin_Eu_Np/_cnt,label='Eu-Npyd')

 #plt.plot(_time,score_ph,label="phenyl")
# plt.plot(_time,score_all,label="sum")
# plt.title("score distribution")
plt.legend()
plt.xlabel("bond (Angs)")
plt.savefig("bond_dist_%s.png"%xtc[:-4],dpi=300)

