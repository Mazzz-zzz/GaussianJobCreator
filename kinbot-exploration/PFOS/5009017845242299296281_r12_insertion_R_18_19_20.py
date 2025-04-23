import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_r12_insertion_R_18_19_20'
logfile = '5009017845242299296281_r12_insertion_R_18_19_20.log'

atom = ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H']
geom = [[0.077309, 0.132052, 0.210786], [-0.501965, -1.224907, -0.260206], [-2.073047, -1.249167, -0.05089], [-2.809932, -1.009184, -1.415494], [-4.358706, -0.979788, -1.178465], [-4.793252, -2.193803, -0.85381], [-4.668808, -0.229611, -0.131183], [-5.331694, -0.402491, -2.679041], [-5.556511, 0.990675, -2.577862], [-4.790011, -1.075206, -3.806621], [-6.695054, -1.123619, -2.352871], [-2.399635, 0.119555, -1.991889], [-2.492725, -1.959382, -2.293768], [-2.43651, -2.422884, 0.474181], [-2.451671, -0.36248, 0.868922], [0.041976, -2.238702, 0.413921], [-0.185362, -1.490985, -1.527237], [1.676355, 0.234885, -0.220928], [1.738273, 1.727204, 0.526128], [3.480849, 0.789707, -0.435388], [3.962359, -0.124005, -0.91453], [3.584798, 1.791875, -1.118416], [4.076257, 1.033862, 0.710335], [2.063066, 1.854874, 1.661198], [1.641843, 2.746884, -0.082459], [2.11613, -0.857927, 0.350602], [1.586483, 0.109779, -1.47573], [-0.596001, 1.146103, -0.301708], [-0.063096, 0.26802, 1.511898], [-6.807853, -2.079514, -2.441352]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_r12_insertion_R_18_19_20', 'label': '5009017845242299296281_r12_insertion_R_18_19_20', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 18 F\n1 28 F\n1 29 F\n2 3 F\n2 16 F\n2 17 F\n3 4 F\n3 14 F\n3 15 F\n4 5 F\n4 12 F\n4 13 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n8 10 F\n8 11 F\n11 30 F\n18 19 F\n18 26 F\n18 27 F\n19 20 F\n19 24 F\n19 25 F\n20 21 F\n20 22 F\n20 23 F\n'}
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
