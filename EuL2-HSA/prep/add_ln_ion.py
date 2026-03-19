import sys
sys.path.append('/Users/hanwei/Desktop/work-new/ff_op')
import PDB
import GRO
import readTop
import topfile
import readFF_new
import connectivity
import geom_opt
import numpy as np


fnm = 'SES_gmx.pdb'
my_pdb = PDB.PDBdata(fnm)

xp = my_pdb.xp
xp1 = xp[my_pdb.index(0,'N1')]
xp2 = xp[my_pdb.index(0,'N2')]
xp3 = xp[my_pdb.index(0,'N3')]

xp_ln = geom_opt.calxyz(xp3,xp2,xp1, 2.65, 56.83, -47.22)


_new_atoms=[]
aid = -1
rid = -1
for _i in my_pdb.atoms:
 _j = [_k for _k in _i]
 _j[0]+=1
 aid =_j[0]
 rid = _j[4]
 _new_atoms.append(_j)


ln_ion=PDB.create_atom(atom_id=aid+1,\
                atom_name = 'EU3',\
                res_name = 'EU3',\
                chain_name = '',\
                res_id=rid+1,\
                x=xp_ln[0],y=xp_ln[1],z=xp_ln[2],\
                seg_name = 'XXXX')

_new_atoms.append(ln_ion)

my_pdb1 = PDB.PDBdata()
my_pdb1.atoms_update(_new_atoms)

my_pdb1.show_pdb(fnm[:-4]+'_c.pdb')

#my_gro = GRO.gro(pdb_obj = my_pdb1, box=[5.0,5.0,5.0])
#my_gro.show(fnm[:-4]+'.gro')



