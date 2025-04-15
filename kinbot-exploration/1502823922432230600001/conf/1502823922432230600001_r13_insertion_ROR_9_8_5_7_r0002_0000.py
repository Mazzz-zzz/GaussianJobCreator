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
            positions=[[np.float64(0.010856), np.float64(2.564688), np.float64(1.910421)], [np.float64(0.173031), np.float64(1.633394), np.float64(2.851845)], [np.float64(-1.240545), np.float64(2.99728), np.float64(1.993043)], [np.float64(0.81984), np.float64(3.575291), np.float64(2.206821)], [np.float64(2.055608), np.float64(1.125765), np.float64(0.065339)], [np.float64(2.818655), np.float64(0.725222), np.float64(1.198311)], [np.float64(0.273463), np.float64(2.098592), np.float64(0.689777)], [np.float64(0.926232), np.float64(0.170394), np.float64(-0.255396)], [np.float64(0.15381), np.float64(0.900282), np.float64(0.231041)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/1502823922432230600001_r13_insertion_ROR_9_8_5_7_r0002_0000', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None'}
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

for c in [[np.int64(7), np.int64(9), np.int64(8), np.int64(5), np.float64(6.418113523169707)]]:
    const.fix_dihedral((c[0]-1, c[1]-1, c[2]-1, c[3]-1), target=c[4])

if os.path.isfile('conf/1502823922432230600001_r13_insertion_ROR_9_8_5_7_r0002_0000_sella.log'):
    os.remove('conf/1502823922432230600001_r13_insertion_ROR_9_8_5_7_r0002_0000_sella.log')

sella_kwargs = {}
opt = Sella(mol, 
            order=0, 
            constraints=const,
            trajectory='conf/1502823922432230600001_r13_insertion_ROR_9_8_5_7_r0002_0000.traj', 
            logfile='conf/1502823922432230600001_r13_insertion_ROR_9_8_5_7_r0002_0000_sella.log',
            **sella_kwargs,
            )

try:
    mol.calc.label = 'conf/1502823922432230600001_r13_insertion_ROR_9_8_5_7_r0002_0000'
    opt.run(fmax=1e-4, steps=100)
    e = mol.get_potential_energy()
    db.write(mol, name='conf/1502823922432230600001_r13_insertion_ROR_9_8_5_7_r0002_0000', 
             data={'energy': e, 'status': 'normal'})
except (RuntimeError, ValueError):
    data = {'status': 'error'}
    db.write(mol, name='conf/1502823922432230600001_r13_insertion_ROR_9_8_5_7_r0002_0000', data=data)
with open('conf/1502823922432230600001_r13_insertion_ROR_9_8_5_7_r0002_0000.log', 'a') as f:
    f.write('done\n')
