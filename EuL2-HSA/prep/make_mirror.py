import sys
sys.path.append('/Users/hanwei/Desktop/work-new/ff_op')
import PDB
import numpy as np
from geom_opt  import calxyz
my_pdb = PDB.PDBdata('SES_gmx_c.pdb')


xp = np.array(my_pdb.xp)

xp[:,2]*=-1.0


_dH = xp[93-1]- xp[89-1]
_dC = xp[92-1]- xp[89-1]
_dH/=np.linalg.norm(_dH)
_dC/=np.linalg.norm(_dC)

xp[92-1] = xp[89-1]+ _dH*1.51
xp[93-1] = xp[89-1]+ _dC*1.0
xp[96-1,:] = np.array(calxyz(xp[92-1],xp[89-1],xp[93-1], 1.0, 109.0 , 180.0))
xp[97-1,:] = np.array(calxyz(xp[92-1],xp[89-1],xp[93-1], 1.0, 109.0 , -60.0))
xp[98-1,:] = np.array(calxyz(xp[92-1],xp[89-1],xp[93-1], 1.0, 109.0 , 60.0))

_dH = xp[82-1]- xp[79-1]
_dC = xp[81-1]- xp[79-1]
_dH/=np.linalg.norm(_dH)
_dC/=np.linalg.norm(_dC)

xp[81-1] = xp[79-1]+ _dH*1.51
xp[82-1] = xp[79-1]+ _dC*1.0
xp[86-1,:] = np.array(calxyz(xp[81-1],xp[79-1],xp[82-1], 1.0, 109.0 , 180.0))
xp[87-1,:] = np.array(calxyz(xp[81-1],xp[79-1],xp[82-1], 1.0, 109.0 , -60.0))
xp[88-1,:] = np.array(calxyz(xp[81-1],xp[79-1],xp[82-1], 1.0, 109.0 , 60.0))


for _i, _j in zip(my_pdb.xp, xp):
 _i[0]= _j[0]
 _i[1]= _j[1]
 _i[2]= _j[2]

my_pdb.update_xp()
my_pdb.show_pdb('SES_gmx_c_m.pdb')

# update_xp()
