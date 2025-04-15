import os
import sys
import shutil

import numpy as np
from ase import Atoms
from ase.db import connect
from sella import Sella, Constraints

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot.stationary_pt import StationaryPoint


db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1502823922432230600001/kinbot.db')
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], 
            positions=[[np.float64(-0.886156), np.float64(-0.389604), np.float64(-0.579583)], [np.float64(-1.675619), np.float64(0.06072), np.float64(-1.54792)], [np.float64(-1.647187), np.float64(-0.966805), np.float64(0.340992)], [np.float64(-0.108907), np.float64(-1.337882), np.float64(-1.105252)], [np.float64(1.297747), np.float64(1.770378), np.float64(-1.042137)], [np.float64(1.686888), np.float64(1.13693), np.float64(-2.255984)], [np.float64(-0.176997), np.float64(0.608544), np.float64(-0.052886)], [np.float64(1.97922), np.float64(1.185044), np.float64(0.176298)], [np.float64(1.010198), np.float64(0.590069), np.float64(0.447675)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/1502823922432230600001_intra_H_migration_7_9_r0001_0000', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None'}
mol.calc = Gaussian(**kwargs)
if 'Gaussian' == 'Gaussian':
    mol.get_potential_energy()
    mol.calc = Gaussian(**kwargs)

const = Constraints(mol)
fix_these = [[idx - 1 for idx in fix] for fix in [[1, 2], [1, 3], [1, 4], [1, 7], [5, 6], [5, 7], [5, 8], [7, 9], [8, 9]]]
for fix in fix_these:
    if len(fix) == 2:
        const.fix_bond(fix)
    elif len(fix) == 4:
        const.fix_dihedral(fix)
    else:
        raise ValueError(f'Unexpected length of fix: {fix}.')

for c in [[np.int64(7), np.int64(9), np.int64(8), np.int64(5), np.float64(0.0)]]:
    const.fix_dihedral((c[0]-1, c[1]-1, c[2]-1, c[3]-1), target=c[4])

if os.path.isfile('conf/1502823922432230600001_intra_H_migration_7_9_r0001_0000_sella.log'):
    os.remove('conf/1502823922432230600001_intra_H_migration_7_9_r0001_0000_sella.log')

sella_kwargs = {}
opt = Sella(mol, 
            order=0, 
            constraints=const,
            trajectory='conf/1502823922432230600001_intra_H_migration_7_9_r0001_0000.traj', 
            logfile='conf/1502823922432230600001_intra_H_migration_7_9_r0001_0000_sella.log',
            **sella_kwargs,
            )

try:
    mol.calc.label = 'conf/1502823922432230600001_intra_H_migration_7_9_r0001_0000'
    opt.run(fmax=1e-4, steps=100)
    e = mol.get_potential_energy()
    db.write(mol, name='conf/1502823922432230600001_intra_H_migration_7_9_r0001_0000', 
             data={'energy': e, 'status': 'normal'})
except (RuntimeError, ValueError):
    data = {'status': 'error'}
    db.write(mol, name='conf/1502823922432230600001_intra_H_migration_7_9_r0001_0000', data=data)
with open('conf/1502823922432230600001_intra_H_migration_7_9_r0001_0000.log', 'a') as f:
    f.write('done\n')
