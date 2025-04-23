import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_r12_insertion_R_26_18_1'
logfile = '5009017845242299296281_r12_insertion_R_26_18_1.log'

atom = ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H']
geom = [[-0.022148, 0.130648, 0.374119], [-0.776557, -0.509049, -0.713266], [-2.345587, -0.335911, -0.436275], [-3.152323, -0.947775, -1.628687], [-4.665857, -1.059273, -1.223818], [-4.839999, -2.093873, -0.413954], [-5.043902, -0.012981, -0.505903], [-5.812411, -1.269806, -2.69484], [-6.255951, 0.007246, -3.11477], [-5.226161, -2.26536, -3.524001], [-7.011197, -1.967374, -1.948298], [-3.005997, -0.203211, -2.719518], [-2.675307, -2.135812, -1.976408], [-2.661527, -0.920928, 0.724034], [-2.647784, 0.9398, -0.255715], [-0.496385, -1.808843, -0.84741], [-0.48622, 0.001383, -1.92819], [1.83903, 0.836992, -0.562364], [2.750542, 1.407289, 0.321454], [4.160607, 0.663863, 0.33214], [4.145052, -0.428249, 1.070716], [4.546407, 0.301433, -0.872295], [5.113523, 1.425892, 0.813195], [2.297552, 1.446759, 1.606208], [2.979522, 2.731555, 0.076067], [1.590464, -0.750728, -0.108305], [1.846937, 0.801514, -1.692929], [-0.340551, 1.301874, 0.682828], [0.19963, -0.483933, 1.421404], [-6.969545, -2.881324, -1.637109]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_r12_insertion_R_26_18_1', 'label': '5009017845242299296281_r12_insertion_R_26_18_1', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 18 F\n1 28 F\n1 29 F\n2 3 F\n2 16 F\n2 17 F\n3 4 F\n3 14 F\n3 15 F\n4 5 F\n4 12 F\n4 13 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n8 10 F\n8 11 F\n11 30 F\n18 19 F\n18 26 F\n18 27 F\n19 20 F\n19 24 F\n19 25 F\n20 21 F\n20 22 F\n20 23 F\n'}
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
