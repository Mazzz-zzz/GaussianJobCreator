import os
import sys
import shutil

import numpy as np
from ase import Atoms
from ase.db import connect
from sella import Sella, Constraints

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot.stationary_pt import StationaryPoint


db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1503305764823040000001/kinbot.db')
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], 
            positions=[[np.float64(-0.060114), np.float64(0.020726), np.float64(-0.043853)], [np.float64(-0.561339), np.float64(0.501262), np.float64(-1.164891)], [np.float64(-0.476131), np.float64(0.72846), np.float64(0.971179)], [np.float64(-0.431402), np.float64(-1.237119), np.float64(0.09953)], [np.float64(1.825269), np.float64(0.042448), np.float64(-0.185863)], [np.float64(2.156797), np.float64(-0.455324), np.float64(1.264867)], [np.float64(3.495277), np.float64(0.164648), np.float64(1.111772)], [np.float64(2.017425), np.float64(1.573076), np.float64(-0.042602)], [np.float64(2.89429), np.float64(1.383872), np.float64(0.591434)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/1503305764823040000001_Intra_RH_Add_Exocyclic_R_7_5_9_r0000_0000', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None'}
mol.calc = Gaussian(**kwargs)
if 'Gaussian' == 'Gaussian':
    mol.get_potential_energy()
    mol.calc = Gaussian(**kwargs)

const = Constraints(mol)
fix_these = [[idx - 1 for idx in fix] for fix in [[1, 2], [1, 3], [1, 4], [1, 5], [5, 6], [5, 7], [5, 8], [6, 7], [7, 9], [8, 9]]]
for fix in fix_these:
    if len(fix) == 2:
        const.fix_bond(fix)
    elif len(fix) == 4:
        const.fix_dihedral(fix)
    else:
        raise ValueError(f'Unexpected length of fix: {fix}.')

for c in [[np.int64(8), np.int64(5), np.int64(6), np.int64(7), np.float64(-45.14265814966824)]]:
    const.fix_dihedral((c[0]-1, c[1]-1, c[2]-1, c[3]-1), target=c[4])

if os.path.isfile('conf/1503305764823040000001_Intra_RH_Add_Exocyclic_R_7_5_9_r0000_0000_sella.log'):
    os.remove('conf/1503305764823040000001_Intra_RH_Add_Exocyclic_R_7_5_9_r0000_0000_sella.log')

sella_kwargs = {}
opt = Sella(mol, 
            order=0, 
            constraints=const,
            trajectory='conf/1503305764823040000001_Intra_RH_Add_Exocyclic_R_7_5_9_r0000_0000.traj', 
            logfile='conf/1503305764823040000001_Intra_RH_Add_Exocyclic_R_7_5_9_r0000_0000_sella.log',
            **sella_kwargs,
            )

try:
    mol.calc.label = 'conf/1503305764823040000001_Intra_RH_Add_Exocyclic_R_7_5_9_r0000_0000'
    opt.run(fmax=1e-4, steps=100)
    e = mol.get_potential_energy()
    db.write(mol, name='conf/1503305764823040000001_Intra_RH_Add_Exocyclic_R_7_5_9_r0000_0000', 
             data={'energy': e, 'status': 'normal'})
except (RuntimeError, ValueError):
    data = {'status': 'error'}
    db.write(mol, name='conf/1503305764823040000001_Intra_RH_Add_Exocyclic_R_7_5_9_r0000_0000', data=data)
with open('conf/1503305764823040000001_Intra_RH_Add_Exocyclic_R_7_5_9_r0000_0000.log', 'a') as f:
    f.write('done\n')
