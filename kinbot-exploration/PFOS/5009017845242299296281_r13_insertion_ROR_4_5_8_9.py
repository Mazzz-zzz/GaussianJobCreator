import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_r13_insertion_ROR_4_5_8_9'
logfile = '5009017845242299296281_r13_insertion_ROR_4_5_8_9.log'

atom = ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H']
geom = [[-2.387047, -1.611784, 1.848054], [-0.87993, -1.599782, 1.430685], [-0.739781, -1.178641, -0.092179], [-0.312961, 0.273047, -0.267005], [1.646772, -0.111461, -0.159633], [1.728376, -1.116993, 0.718409], [1.744273, -0.719652, -1.342778], [2.301858, 1.163278, 0.060397], [0.57261, 2.062217, -0.388567], [2.985067, 1.396518, 1.200245], [3.263389, 1.51771, -1.068032], [-0.54342, 0.630443, -1.444156], [-0.863937, 0.962978, 0.612585], [0.037294, -2.088372, -0.71962], [-1.905458, -1.354595, -0.724769], [-0.339971, -2.803758, 1.608376], [-0.170301, -0.812419, 2.228352], [-2.488545, -1.79019, 3.406591], [-3.982752, -2.080548, 3.790883], [-4.143146, -1.991732, 5.346791], [-3.25592, -2.729197, 5.996369], [-4.01269, -0.759594, 5.800803], [-5.329858, -2.410061, 5.742398], [-4.363768, -3.274167, 3.350992], [-4.814127, -1.230607, 3.19194], [-1.698629, -2.767014, 3.833303], [-2.0368, -0.704952, 4.030963], [-3.004356, -0.490368, 1.470448], [-3.045937, -2.577791, 1.212746], [4.223566, 1.53832, -0.974153]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_r13_insertion_ROR_4_5_8_9', 'label': '5009017845242299296281_r13_insertion_ROR_4_5_8_9', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 18 F\n1 28 F\n1 29 F\n2 3 F\n2 16 F\n2 17 F\n3 4 F\n3 14 F\n3 15 F\n4 5 F\n4 12 F\n4 13 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n8 10 F\n8 11 F\n11 30 F\n18 19 F\n18 26 F\n18 27 F\n19 20 F\n19 24 F\n19 25 F\n20 21 F\n20 22 F\n20 23 F\n'}
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
