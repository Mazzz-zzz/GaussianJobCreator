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
            positions=[[np.float64(3.460359), np.float64(0.65917), np.float64(0.127396)], [np.float64(3.886376), np.float64(0.113531), np.float64(-0.994967)], [np.float64(3.46589), np.float64(-0.227426), np.float64(1.085579)], [np.float64(4.240914), np.float64(1.674), np.float64(0.446081)], [np.float64(1.72661), np.float64(1.349909), np.float64(-0.176353)], [np.float64(1.424726), np.float64(1.826576), np.float64(1.287504)], [np.float64(-0.017657), np.float64(1.763363), np.float64(0.947564)], [np.float64(0.978814), np.float64(-0.006277), np.float64(-0.224051)], [np.float64(0.159101), np.float64(0.44889), np.float64(0.348709)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/1503305764823040000001_r13_insertion_ROR_9_8_5_7_r0000_0000', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None'}
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

for c in [[np.int64(6), np.int64(7), np.int64(9), np.int64(8), np.float64(-45.271506705625576)]]:
    const.fix_dihedral((c[0]-1, c[1]-1, c[2]-1, c[3]-1), target=c[4])

if os.path.isfile('conf/1503305764823040000001_r13_insertion_ROR_9_8_5_7_r0000_0000_sella.log'):
    os.remove('conf/1503305764823040000001_r13_insertion_ROR_9_8_5_7_r0000_0000_sella.log')

sella_kwargs = {}
opt = Sella(mol, 
            order=0, 
            constraints=const,
            trajectory='conf/1503305764823040000001_r13_insertion_ROR_9_8_5_7_r0000_0000.traj', 
            logfile='conf/1503305764823040000001_r13_insertion_ROR_9_8_5_7_r0000_0000_sella.log',
            **sella_kwargs,
            )

try:
    mol.calc.label = 'conf/1503305764823040000001_r13_insertion_ROR_9_8_5_7_r0000_0000'
    opt.run(fmax=1e-4, steps=100)
    e = mol.get_potential_energy()
    db.write(mol, name='conf/1503305764823040000001_r13_insertion_ROR_9_8_5_7_r0000_0000', 
             data={'energy': e, 'status': 'normal'})
except (RuntimeError, ValueError):
    data = {'status': 'error'}
    db.write(mol, name='conf/1503305764823040000001_r13_insertion_ROR_9_8_5_7_r0000_0000', data=data)
with open('conf/1503305764823040000001_r13_insertion_ROR_9_8_5_7_r0000_0000.log', 'a') as f:
    f.write('done\n')
