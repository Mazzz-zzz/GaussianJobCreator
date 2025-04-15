import os
import sys
import shutil

import numpy as np
from ase import Atoms
from ase.db import connect
from sella import Sella, Constraints

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot.stationary_pt import StationaryPoint


db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1502984102222700000001/kinbot.db')
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], 
            positions=[[np.float64(1.720925), np.float64(2.171827), np.float64(1.318126)], [np.float64(0.977157), np.float64(2.141359), np.float64(2.386677)], [np.float64(0.530488), np.float64(2.450857), np.float64(0.202863)], [np.float64(2.358916), np.float64(3.294641), np.float64(1.249307)], [np.float64(1.159856), np.float64(-0.103292), np.float64(-0.078883)], [np.float64(0.616077), np.float64(-0.972427), np.float64(0.904015)], [np.float64(2.309195), np.float64(1.130308), np.float64(0.941148)], [np.float64(2.067315), np.float64(-0.552413), np.float64(-1.072038)], [np.float64(0.490277), np.float64(1.282511), np.float64(-0.22401)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/1502984102222700000001_r12_insertion_R_7_1_3_r0000_0000', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None'}
mol.calc = Gaussian(**kwargs)
if 'Gaussian' == 'Gaussian':
    mol.get_potential_energy()
    mol.calc = Gaussian(**kwargs)

const = Constraints(mol)
fix_these = [[idx - 1 for idx in fix] for fix in [[1, 2], [1, 3], [1, 4], [1, 7], [3, 9], [5, 6], [5, 7], [5, 8], [5, 9]]]
for fix in fix_these:
    if len(fix) == 2:
        const.fix_bond(fix)
    elif len(fix) == 4:
        const.fix_dihedral(fix)
    else:
        raise ValueError(f'Unexpected length of fix: {fix}.')

for c in [[np.int64(9), np.int64(5), np.int64(7), np.int64(1), np.float64(-20.447545334510995)]]:
    const.fix_dihedral((c[0]-1, c[1]-1, c[2]-1, c[3]-1), target=c[4])

if os.path.isfile('conf/1502984102222700000001_r12_insertion_R_7_1_3_r0000_0000_sella.log'):
    os.remove('conf/1502984102222700000001_r12_insertion_R_7_1_3_r0000_0000_sella.log')

sella_kwargs = {}
opt = Sella(mol, 
            order=0, 
            constraints=const,
            trajectory='conf/1502984102222700000001_r12_insertion_R_7_1_3_r0000_0000.traj', 
            logfile='conf/1502984102222700000001_r12_insertion_R_7_1_3_r0000_0000_sella.log',
            **sella_kwargs,
            )

try:
    mol.calc.label = 'conf/1502984102222700000001_r12_insertion_R_7_1_3_r0000_0000'
    opt.run(fmax=1e-4, steps=100)
    e = mol.get_potential_energy()
    db.write(mol, name='conf/1502984102222700000001_r12_insertion_R_7_1_3_r0000_0000', 
             data={'energy': e, 'status': 'normal'})
except (RuntimeError, ValueError):
    data = {'status': 'error'}
    db.write(mol, name='conf/1502984102222700000001_r12_insertion_R_7_1_3_r0000_0000', data=data)
with open('conf/1502984102222700000001_r12_insertion_R_7_1_3_r0000_0000.log', 'a') as f:
    f.write('done\n')
