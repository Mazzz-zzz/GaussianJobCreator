import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_r13_insertion_RSR_12_4_5_8'
logfile = '5009017845242299296281_r13_insertion_RSR_12_4_5_8.log'

atom = ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H']
geom = [[-0.022325, -2.440179, -2.286253], [0.59612, -1.93869, -0.943608], [1.427937, -0.615808, -1.223295], [1.858636, 0.02082, 0.078769], [2.381741, 1.373154, 0.086158], [3.43743, 1.546187, 0.87568], [2.787872, 1.817402, -1.093347], [0.946237, 2.575454, 0.788872], [0.051026, 2.961779, -0.206476], [0.867193, 2.341302, 2.168727], [1.974715, 3.751254, 0.754591], [0.089193, 0.775775, 0.62568], [1.997763, -0.724143, 1.084354], [2.511678, -0.913071, -1.971038], [0.734989, 0.230522, -1.99507], [1.387261, -2.85461, -0.39355], [-0.345984, -1.726948, -0.029518], [-0.605685, -3.884837, -2.0689], [-1.462668, -4.284055, -3.318738], [-1.825708, -5.805539, -3.231871], [-0.772693, -6.584587, -3.401484], [-2.359042, -6.139273, -2.074113], [-2.697045, -6.15328, -4.157133], [-0.82111, -4.024869, -4.452541], [-2.567838, -3.545206, -3.392083], [0.377369, -4.75632, -1.859737], [-1.338395, -3.945398, -0.95808], [-0.961703, -1.601814, -2.715705], [0.882959, -2.447181, -3.262071], [2.687877, 3.882768, 1.39594]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_r13_insertion_RSR_12_4_5_8', 'label': '5009017845242299296281_r13_insertion_RSR_12_4_5_8', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 18 F\n1 28 F\n1 29 F\n2 3 F\n2 16 F\n2 17 F\n3 4 F\n3 14 F\n3 15 F\n4 5 F\n4 12 F\n4 13 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n8 10 F\n8 11 F\n11 30 F\n18 19 F\n18 26 F\n18 27 F\n19 20 F\n19 24 F\n19 25 F\n20 21 F\n20 22 F\n20 23 F\n'}
Gaussian.command = 'g16 < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy() # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    freq = reader_gauss.read_freq(logfile, ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'])
    zpe = reader_gauss.read_zpe(logfile)
    db.write(mol, name=label, data={'energy': e,'frequencies': np.asarray(freq), 'zpe':zpe, 'status': 'normal'})
except RuntimeError:
    try:
        iowait(logfile, 'gauss')
        mol.positions = reader_gauss.read_geom(logfile, mol)
        kwargs = reader_gauss.correct_kwargs(logfile, kwargs)
        mol.calc = Gaussian(**kwargs)
        e = mol.get_potential_energy()  # use the Gaussian optimizer
        iowait(logfile, 'gauss')
        mol.positions = reader_gauss.read_geom(logfile, mol)
        freq = reader_gauss.read_freq(logfile, ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'])
        zpe = reader_gauss.read_zpe(logfile)
        db.write(mol, name=label, data={'energy': e,
                                         'frequencies': np.asarray(freq),
                                         'zpe': zpe, 'status': 'normal'})
    except RuntimeError:
        db.write(mol, name=label, data={'status': 'error'})

with open(logfile,'a') as f:
    f.write('done\n')
