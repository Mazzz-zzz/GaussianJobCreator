import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_r12_insertion_R_12_4_3'
logfile = '5009017845242299296281_r12_insertion_R_12_4_3.log'

atom = ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H']
geom = [[0.205197, 0.470847, 0.129328], [-0.771857, -0.382106, -0.773799], [-2.176282, -0.02471, -0.448744], [-3.364102, -1.756216, -1.105197], [-4.747076, -1.617691, -0.995747], [-5.398264, -2.771815, -0.67562], [-5.123752, -0.801514, 0.024589], [-5.578249, -0.982255, -2.597343], [-5.60737, 0.427068, -2.618812], [-5.118336, -1.844362, -3.62785], [-7.044871, -1.454952, -2.280865], [-2.811002, -0.578333, -2.151883], [-2.72831, -2.669086, -1.341592], [-2.515486, -0.169211, 0.753156], [-2.658049, 1.030254, -0.85336], [-0.518851, -1.694655, -0.583413], [-0.470473, -0.196852, -2.074312], [1.674191, 0.300916, -0.389863], [2.66168, 0.905997, 0.668735], [4.09436, 0.990818, 0.05284], [4.494903, -0.153648, -0.46977], [4.181967, 1.892344, -0.902869], [4.996296, 1.318523, 0.956269], [2.669746, 0.183064, 1.779382], [2.258336, 2.115571, 1.07054], [1.961061, -0.974236, -0.621926], [1.821626, 0.893053, -1.574939], [-0.149796, 1.754031, 0.129542], [0.094281, 0.102801, 1.405123], [-7.297669, -2.391825, -2.26349]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_r12_insertion_R_12_4_3', 'label': '5009017845242299296281_r12_insertion_R_12_4_3', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 18 F\n1 28 F\n1 29 F\n2 3 F\n2 16 F\n2 17 F\n3 4 F\n3 14 F\n3 15 F\n4 5 F\n4 12 F\n4 13 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n8 10 F\n8 11 F\n11 30 F\n18 19 F\n18 26 F\n18 27 F\n19 20 F\n19 24 F\n19 25 F\n20 21 F\n20 22 F\n20 23 F\n'}
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
