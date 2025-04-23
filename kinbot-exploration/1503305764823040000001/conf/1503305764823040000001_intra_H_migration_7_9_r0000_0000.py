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
            positions=[[np.float64(2.945984), np.float64(0.857347), np.float64(-1.619794)], [np.float64(3.937172), np.float64(1.702032), np.float64(-1.412492)], [np.float64(2.1516), np.float64(1.334443), np.float64(-2.558878)], [np.float64(3.417783), np.float64(-0.305015), np.float64(-1.98101)], [np.float64(1.934042), np.float64(0.738096), np.float64(-0.026878)], [np.float64(0.922503), np.float64(-0.366052), np.float64(-0.495083)], [np.float64(0.682187), np.float64(-0.67274), np.float64(0.935969)], [np.float64(2.950928), np.float64(-0.096505), np.float64(0.791813)], [np.float64(2.114912), np.float64(-0.682369), np.float64(1.196549)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/1503305764823040000001_intra_H_migration_7_9_r0000_0000', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None'}
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

for c in [[np.int64(5), np.int64(6), np.int64(7), np.int64(9), np.float64(-45.17574531437602)]]:
    const.fix_dihedral((c[0]-1, c[1]-1, c[2]-1, c[3]-1), target=c[4])

if os.path.isfile('conf/1503305764823040000001_intra_H_migration_7_9_r0000_0000_sella.log'):
    os.remove('conf/1503305764823040000001_intra_H_migration_7_9_r0000_0000_sella.log')

sella_kwargs = {}
opt = Sella(mol, 
            order=0, 
            constraints=const,
            trajectory='conf/1503305764823040000001_intra_H_migration_7_9_r0000_0000.traj', 
            logfile='conf/1503305764823040000001_intra_H_migration_7_9_r0000_0000_sella.log',
            **sella_kwargs,
            )

try:
    mol.calc.label = 'conf/1503305764823040000001_intra_H_migration_7_9_r0000_0000'
    opt.run(fmax=1e-4, steps=100)
    e = mol.get_potential_energy()
    db.write(mol, name='conf/1503305764823040000001_intra_H_migration_7_9_r0000_0000', 
             data={'energy': e, 'status': 'normal'})
except (RuntimeError, ValueError):
    data = {'status': 'error'}
    db.write(mol, name='conf/1503305764823040000001_intra_H_migration_7_9_r0000_0000', data=data)
with open('conf/1503305764823040000001_intra_H_migration_7_9_r0000_0000.log', 'a') as f:
    f.write('done\n')
