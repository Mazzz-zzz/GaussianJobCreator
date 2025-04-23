import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_intra_H_migration_8_30'
logfile = '5009017845242299296281_intra_H_migration_8_30.log'

atom = ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H']
geom = [[0.012454, 0.082174, 0.138963], [-0.647901, -1.289356, -0.214059], [-2.207482, -1.211627, 0.010425], [-2.955207, -1.096991, -1.366691], [-4.500343, -1.034964, -1.108468], [-4.940059, -2.228248, -0.723549], [-4.785835, -0.234646, -0.093661], [-5.470164, -0.482324, -2.618544], [-5.451676, 0.931508, -2.665167], [-4.930837, -1.22832, -3.703762], [-6.904657, -0.993179, -2.220098], [-2.537825, -0.034024, -2.047523], [-2.660606, -2.127122, -2.150465], [-2.615259, -2.299536, 0.671181], [-2.530397, -0.209697, 0.819273], [-0.138593, -2.263596, 0.536547], [-0.34645, -1.660682, -1.455933], [1.532093, 0.037476, -0.263392], [2.257087, 1.291823, 0.336645], [3.690561, 1.406279, -0.281578], [4.379866, 0.283413, -0.186603], [3.667418, 1.727205, -1.558841], [4.410668, 2.332428, 0.321126], [2.325792, 1.212545, 1.660817], [1.566707, 2.408127, 0.112161], [2.110774, -1.08251, 0.158912], [1.662907, -0.005098, -1.587399], [-0.600005, 1.088731, -0.479612], [-0.130041, 0.356158, 1.433608], [-6.229136, -1.754128, -2.856196]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_intra_H_migration_8_30', 'label': '5009017845242299296281_intra_H_migration_8_30', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': ''}
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
