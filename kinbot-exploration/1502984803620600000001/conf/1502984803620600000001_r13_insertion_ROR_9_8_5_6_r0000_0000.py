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
            positions=[[np.float64(3.414126), np.float64(0.594216), np.float64(1.245456)], [np.float64(3.943004), np.float64(-0.447201), np.float64(0.667839)], [np.float64(2.967882), np.float64(0.302882), np.float64(2.437884)], [np.float64(4.303163), np.float64(1.554965), np.float64(1.321116)], [np.float64(1.929722), np.float64(1.397657), np.float64(0.169322)], [np.float64(-0.091791), np.float64(1.388931), np.float64(-0.297873)], [np.float64(2.516987), np.float64(1.711118), np.float64(-1.100645)], [np.float64(1.261482), np.float64(0.059468), np.float64(0.304379)], [np.float64(-0.573445), np.float64(1.34165), np.float64(0.540934)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/1502984803620600000001_r13_insertion_ROR_9_8_5_6_r0000_0000', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None'}
mol.calc = Gaussian(**kwargs)
if 'Gaussian' == 'Gaussian':
    mol.get_potential_energy()
    mol.calc = Gaussian(**kwargs)

const = Constraints(mol)
fix_these = [[idx - 1 for idx in fix] for fix in [[1, 2], [1, 3], [1, 4], [1, 5], [5, 6], [5, 7], [5, 8], [6, 8], [6, 9], [8, 9]]]
for fix in fix_these:
    if len(fix) == 2:
        const.fix_bond(fix)
    elif len(fix) == 4:
        const.fix_dihedral(fix)
    else:
        raise ValueError(f'Unexpected length of fix: {fix}.')

for c in [[np.int64(5), np.int64(6), np.int64(9), np.int64(8), np.float64(-50.714414077188735)]]:
    const.fix_dihedral((c[0]-1, c[1]-1, c[2]-1, c[3]-1), target=c[4])

if os.path.isfile('conf/1502984803620600000001_r13_insertion_ROR_9_8_5_6_r0000_0000_sella.log'):
    os.remove('conf/1502984803620600000001_r13_insertion_ROR_9_8_5_6_r0000_0000_sella.log')

sella_kwargs = {}
opt = Sella(mol, 
            order=0, 
            constraints=const,
            trajectory='conf/1502984803620600000001_r13_insertion_ROR_9_8_5_6_r0000_0000.traj', 
            logfile='conf/1502984803620600000001_r13_insertion_ROR_9_8_5_6_r0000_0000_sella.log',
            **sella_kwargs,
            )

try:
    mol.calc.label = 'conf/1502984803620600000001_r13_insertion_ROR_9_8_5_6_r0000_0000'
    opt.run(fmax=1e-4, steps=100)
    e = mol.get_potential_energy()
    db.write(mol, name='conf/1502984803620600000001_r13_insertion_ROR_9_8_5_6_r0000_0000', 
             data={'energy': e, 'status': 'normal'})
except (RuntimeError, ValueError):
    data = {'status': 'error'}
    db.write(mol, name='conf/1502984803620600000001_r13_insertion_ROR_9_8_5_6_r0000_0000', data=data)
with open('conf/1502984803620600000001_r13_insertion_ROR_9_8_5_6_r0000_0000.log', 'a') as f:
    f.write('done\n')
