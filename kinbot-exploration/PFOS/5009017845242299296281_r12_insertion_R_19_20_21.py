import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_r12_insertion_R_19_20_21'
logfile = '5009017845242299296281_r12_insertion_R_19_20_21.log'

atom = ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H']
geom = [[0.094172, -0.063735, 0.317651], [-0.595177, -1.394893, -0.104296], [-2.166361, -1.253714, 0.022285], [-2.807195, -1.057267, -1.396948], [-4.3559, -0.857864, -1.24694], [-4.921849, -1.9833, -0.83592], [-4.632872, 0.018384, -0.29461], [-5.194394, -0.33839, -2.845774], [-5.269633, 1.072223, -2.891281], [-4.672141, -1.175183, -3.874265], [-6.637043, -0.87997, -2.531832], [-2.258537, -0.026931, -2.031609], [-2.546221, -2.100519, -2.176659], [-2.664859, -2.339056, 0.612544], [-2.498136, -0.263982, 0.846469], [-0.180947, -2.404208, 0.661273], [-0.239607, -1.762984, -1.329239], [1.653602, -0.164192, 0.015238], [2.388662, 1.106664, 0.491927], [2.93679, 2.074179, -0.754003], [3.940288, 0.213634, -0.14474], [2.638955, 1.894546, -1.91079], [3.829348, 2.87368, -0.570289], [3.1456, 1.027965, 1.550776], [1.46893, 1.978556, 0.988869], [2.149776, -1.242494, 0.599934], [1.833948, -0.373935, -1.286709], [-0.433138, 0.978662, -0.317842], [-0.117976, 0.188007, 1.609138], [-6.854235, -1.822071, -2.554231]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_r12_insertion_R_19_20_21', 'label': '5009017845242299296281_r12_insertion_R_19_20_21', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 18 F\n1 28 F\n1 29 F\n2 3 F\n2 16 F\n2 17 F\n3 4 F\n3 14 F\n3 15 F\n4 5 F\n4 12 F\n4 13 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n8 10 F\n8 11 F\n11 30 F\n18 19 F\n18 26 F\n18 27 F\n19 20 F\n19 24 F\n19 25 F\n20 21 F\n20 22 F\n20 23 F\n'}
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
