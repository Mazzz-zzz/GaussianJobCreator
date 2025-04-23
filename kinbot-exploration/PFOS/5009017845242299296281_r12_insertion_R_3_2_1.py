import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_r12_insertion_R_3_2_1'
logfile = '5009017845242299296281_r12_insertion_R_3_2_1.log'

atom = ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H']
geom = [[-0.385401, -0.227342, 0.648449], [-0.422976, -2.375936, 0.177217], [-1.800134, -1.435045, 0.261095], [-2.278879, -0.78571, -1.068393], [-3.855485, -0.811756, -1.126864], [-4.282717, -2.039675, -1.409911], [-4.385092, -0.542729, 0.057076], [-4.568432, 0.346673, -2.4237], [-4.83041, 1.607525, -1.838037], [-3.840898, 0.117187, -3.623691], [-5.95489, -0.372407, -2.620285], [-1.863774, 0.443189, -1.305303], [-1.787501, -1.482357, -2.10856], [-2.516197, -2.646849, 0.229822], [-2.357599, -1.038419, 1.347446], [0.088647, -2.843497, 1.13499], [0.07441, -2.674642, -0.84801], [0.778597, -0.071608, -0.319259], [1.414855, 1.37242, -0.189145], [2.784944, 1.407446, -0.943614], [3.722739, 0.714648, -0.321769], [2.715139, 0.923428, -2.166625], [3.253911, 2.632837, -1.060224], [1.57754, 1.712311, 1.079792], [0.592742, 2.294882, -0.67537], [1.733312, -0.99351, -0.101291], [0.44977, -0.291755, -1.59366], [-1.002283, 0.758879, 0.762731], [0.0188, -0.563215, 1.757504], [-6.038887, -1.213153, -3.094296]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_r12_insertion_R_3_2_1', 'label': '5009017845242299296281_r12_insertion_R_3_2_1', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 18 F\n1 28 F\n1 29 F\n2 3 F\n2 16 F\n2 17 F\n3 4 F\n3 14 F\n3 15 F\n4 5 F\n4 12 F\n4 13 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n8 10 F\n8 11 F\n11 30 F\n18 19 F\n18 26 F\n18 27 F\n19 20 F\n19 24 F\n19 25 F\n20 21 F\n20 22 F\n20 23 F\n'}
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
