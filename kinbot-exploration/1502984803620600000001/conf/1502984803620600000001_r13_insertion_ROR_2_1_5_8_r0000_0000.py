import os
import sys
import shutil

import numpy as np
from ase import Atoms
from ase.db import connect
from sella import Sella, Constraints

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot.stationary_pt import StationaryPoint


db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1502984803620600000001/kinbot.db')
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], 
            positions=[[np.float64(1.545794), np.float64(0.107648), np.float64(-0.085506)], [np.float64(-0.402699), np.float64(0.501464), np.float64(0.185534)], [np.float64(1.688131), np.float64(-0.65989), np.float64(0.899926)], [np.float64(1.456075), np.float64(-0.513791), np.float64(-1.171988)], [np.float64(2.056436), np.float64(1.986146), np.float64(0.00044)], [np.float64(3.244694), np.float64(1.860497), np.float64(0.775864)], [np.float64(2.141694), np.float64(2.288084), np.float64(-1.387002)], [np.float64(0.880682), np.float64(2.479777), np.float64(0.701647)], [np.float64(-0.184033), np.float64(1.419948), np.float64(0.452889)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/1502984803620600000001_r13_insertion_ROR_2_1_5_8_r0000_0000', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None'}
mol.calc = Gaussian(**kwargs)
if 'Gaussian' == 'Gaussian':
    mol.get_potential_energy()
    mol.calc = Gaussian(**kwargs)

const = Constraints(mol)
fix_these = [[idx - 1 for idx in fix] for fix in [[1, 2], [1, 3], [1, 4], [1, 5], [2, 9], [5, 6], [5, 7], [5, 8], [8, 9]]]
for fix in fix_these:
    if len(fix) == 2:
        const.fix_bond(fix)
    elif len(fix) == 4:
        const.fix_dihedral(fix)
    else:
        raise ValueError(f'Unexpected length of fix: {fix}.')

for c in [[np.int64(2), np.int64(9), np.int64(8), np.int64(5), np.float64(-18.228869487284417)]]:
    const.fix_dihedral((c[0]-1, c[1]-1, c[2]-1, c[3]-1), target=c[4])

if os.path.isfile('conf/1502984803620600000001_r13_insertion_ROR_2_1_5_8_r0000_0000_sella.log'):
    os.remove('conf/1502984803620600000001_r13_insertion_ROR_2_1_5_8_r0000_0000_sella.log')

sella_kwargs = {}
opt = Sella(mol, 
            order=0, 
            constraints=const,
            trajectory='conf/1502984803620600000001_r13_insertion_ROR_2_1_5_8_r0000_0000.traj', 
            logfile='conf/1502984803620600000001_r13_insertion_ROR_2_1_5_8_r0000_0000_sella.log',
            **sella_kwargs,
            )

try:
    mol.calc.label = 'conf/1502984803620600000001_r13_insertion_ROR_2_1_5_8_r0000_0000'
    opt.run(fmax=1e-4, steps=100)
    e = mol.get_potential_energy()
    db.write(mol, name='conf/1502984803620600000001_r13_insertion_ROR_2_1_5_8_r0000_0000', 
             data={'energy': e, 'status': 'normal'})
except (RuntimeError, ValueError):
    data = {'status': 'error'}
    db.write(mol, name='conf/1502984803620600000001_r13_insertion_ROR_2_1_5_8_r0000_0000', data=data)
with open('conf/1502984803620600000001_r13_insertion_ROR_2_1_5_8_r0000_0000.log', 'a') as f:
    f.write('done\n')
