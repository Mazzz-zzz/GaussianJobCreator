import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_r13_insertion_ROR_6_5_8_9'
logfile = '5009017845242299296281_r13_insertion_ROR_6_5_8_9.log'

atom = ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H']
geom = [[2.411196, -4.308818, -1.464847], [1.503641, -3.234557, -2.144002], [1.032274, -2.172741, -1.073625], [1.860564, -0.849444, -1.228251], [1.360131, 0.202927, -0.162412], [-0.216388, 0.851124, -1.208561], [0.639153, -0.380212, 0.802131], [1.933245, 1.52587, -0.00792], [0.672678, 2.638768, -1.090707], [2.323483, 1.884425, 1.239291], [3.228961, 1.568492, -0.819318], [3.157103, -1.079608, -1.070314], [1.757987, -0.357257, -2.455426], [-0.272827, -1.935806, -1.226137], [1.127464, -2.65567, 0.158098], [0.441175, -3.801778, -2.702898], [2.128926, -2.654737, -3.167989], [3.055103, -5.214615, -2.577726], [3.734039, -6.454386, -1.898091], [4.601434, -7.216177, -2.955524], [3.931015, -7.499182, -4.05819], [5.66404, -6.532463, -3.324569], [5.04351, -8.367357, -2.486405], [2.822357, -7.265791, -1.369431], [4.485899, -6.088792, -0.859651], [2.146093, -5.606546, -3.465961], [3.935339, -4.517799, -3.293496], [3.360076, -3.743454, -0.723357], [1.711352, -5.038792, -0.60049], [4.109857, 1.549392, -0.422321]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_r13_insertion_ROR_6_5_8_9', 'label': '5009017845242299296281_r13_insertion_ROR_6_5_8_9', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 18 F\n1 28 F\n1 29 F\n2 3 F\n2 16 F\n2 17 F\n3 4 F\n3 14 F\n3 15 F\n4 5 F\n4 12 F\n4 13 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n8 10 F\n8 11 F\n11 30 F\n18 19 F\n18 26 F\n18 27 F\n19 20 F\n19 24 F\n19 25 F\n20 21 F\n20 22 F\n20 23 F\n'}
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
