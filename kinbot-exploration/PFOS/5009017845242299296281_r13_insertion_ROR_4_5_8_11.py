import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_r13_insertion_ROR_4_5_8_11'
logfile = '5009017845242299296281_r13_insertion_ROR_4_5_8_11.log'

atom = ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H']
geom = [[-2.639475, -2.473352, 1.19191], [-2.163336, -1.093972, 0.641261], [-0.58125, -0.994153, 0.778038], [-0.145785, 0.35535, 0.215775], [1.810853, -0.058605, 0.228611], [1.920285, -0.799148, 1.342522], [1.797281, -0.970457, -0.752807], [2.566336, 1.169961, 0.078965], [2.98476, 1.833468, 1.289605], [3.289923, 1.384511, -1.160327], [0.811307, 1.998329, -0.404481], [-0.509399, 0.485885, -0.978771], [-0.483088, 1.266322, 1.00217], [-0.246704, -1.129309, 2.073013], [-0.021801, -2.047374, 0.182817], [-2.747423, -0.083778, 1.281791], [-2.534183, -0.933374, -0.627876], [-4.213457, -2.496761, 1.214646], [-4.69943, -3.949247, 1.539653], [-6.246859, -3.930304, 1.801374], [-6.560491, -3.305422, 2.920332], [-6.921315, -3.341682, 0.836558], [-6.737004, -5.15087, 1.907583], [-4.060386, -4.456169, 2.586566], [-4.408363, -4.782339, 0.543728], [-4.680787, -1.621702, 2.100856], [-4.716391, -2.083355, 0.050648], [-2.158957, -3.469777, 0.451414], [-2.155393, -2.705119, 2.408332], [0.857009, 2.631389, -1.138574]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_r13_insertion_ROR_4_5_8_11', 'label': '5009017845242299296281_r13_insertion_ROR_4_5_8_11', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 18 F\n1 28 F\n1 29 F\n2 3 F\n2 16 F\n2 17 F\n3 4 F\n3 14 F\n3 15 F\n4 5 F\n4 12 F\n4 13 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n8 10 F\n8 11 F\n11 30 F\n18 19 F\n18 26 F\n18 27 F\n19 20 F\n19 24 F\n19 25 F\n20 21 F\n20 22 F\n20 23 F\n'}
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
