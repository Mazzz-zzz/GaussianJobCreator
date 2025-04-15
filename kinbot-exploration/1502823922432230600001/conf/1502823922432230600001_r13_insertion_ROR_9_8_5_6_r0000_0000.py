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
            positions=[[np.float64(2.08039), np.float64(1.675295), np.float64(-2.805191)], [np.float64(3.055592), np.float64(1.415208), np.float64(-3.64292)], [np.float64(2.458923), np.float64(2.618574), np.float64(-1.975331)], [np.float64(1.029664), np.float64(2.106845), np.float64(-3.483536)], [np.float64(0.436244), np.float64(0.301303), np.float64(-1.257538)], [np.float64(1.287027), np.float64(0.795789), np.float64(0.951023)], [np.float64(1.794802), np.float64(0.507461), np.float64(-2.141117)], [np.float64(0.144451), np.float64(1.514008), np.float64(-0.523916)], [np.float64(1.660645), np.float64(1.181094), np.float64(1.757558)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/1502823922432230600001_r13_insertion_ROR_9_8_5_6_r0000_0000', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None'}
mol.calc = Gaussian(**kwargs)
if 'Gaussian' == 'Gaussian':
    mol.get_potential_energy()
    mol.calc = Gaussian(**kwargs)

const = Constraints(mol)
fix_these = [[idx - 1 for idx in fix] for fix in [[1, 2], [1, 3], [1, 4], [1, 7], [5, 6], [5, 7], [5, 8], [6, 8], [6, 9], [8, 9]]]
for fix in fix_these:
    if len(fix) == 2:
        const.fix_bond(fix)
    elif len(fix) == 4:
        const.fix_dihedral(fix)
    else:
        raise ValueError(f'Unexpected length of fix: {fix}.')

for c in [[np.int64(5), np.int64(6), np.int64(9), np.int64(8), np.float64(-29.926085026017745)]]:
    const.fix_dihedral((c[0]-1, c[1]-1, c[2]-1, c[3]-1), target=c[4])

if os.path.isfile('conf/1502823922432230600001_r13_insertion_ROR_9_8_5_6_r0000_0000_sella.log'):
    os.remove('conf/1502823922432230600001_r13_insertion_ROR_9_8_5_6_r0000_0000_sella.log')

sella_kwargs = {}
opt = Sella(mol, 
            order=0, 
            constraints=const,
            trajectory='conf/1502823922432230600001_r13_insertion_ROR_9_8_5_6_r0000_0000.traj', 
            logfile='conf/1502823922432230600001_r13_insertion_ROR_9_8_5_6_r0000_0000_sella.log',
            **sella_kwargs,
            )

try:
    mol.calc.label = 'conf/1502823922432230600001_r13_insertion_ROR_9_8_5_6_r0000_0000'
    opt.run(fmax=1e-4, steps=100)
    e = mol.get_potential_energy()
    db.write(mol, name='conf/1502823922432230600001_r13_insertion_ROR_9_8_5_6_r0000_0000', 
             data={'energy': e, 'status': 'normal'})
except (RuntimeError, ValueError):
    data = {'status': 'error'}
    db.write(mol, name='conf/1502823922432230600001_r13_insertion_ROR_9_8_5_6_r0000_0000', data=data)
with open('conf/1502823922432230600001_r13_insertion_ROR_9_8_5_6_r0000_0000.log', 'a') as f:
    f.write('done\n')
