import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_Intra_Diels_alder_R_9_10'
logfile = '5009017845242299296281_Intra_Diels_alder_R_9_10.log'

atom = ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H']
geom = [[0.023739, 0.020744, 0.175969], [-0.619318, -1.335716, -0.259419], [-2.181686, -1.285423, -0.052307], [-2.903342, -1.025495, -1.421318], [-4.449741, -0.889766, -1.164151], [-4.918205, -2.02888, -0.647002], [-4.68798, 0.014623, -0.219033], [-5.4205, -0.503311, -2.665346], [-4.859303, -0.808546, -3.863383], [-5.468666, 0.750931, -3.201717], [-6.72709, -1.152258, -2.226946], [-2.423546, 0.060096, -2.016912], [-2.661342, -2.015663, -2.270312], [-2.596553, -2.44082, 0.479401], [-2.520741, -0.381009, 0.860409], [-0.109291, -2.341898, 0.447029], [-0.303245, -1.642317, -1.515201], [1.530921, 0.041187, -0.275327], [2.249711, 1.259715, 0.399905], [3.665531, 1.446212, -0.240992], [4.374607, 0.331792, -0.247589], [3.606001, 1.862721, -1.488943], [4.385083, 2.336392, 0.414389], [2.352028, 1.083775, 1.712584], [1.535434, 2.376856, 0.275746], [2.142801, -1.096301, 0.039792], [1.620042, 0.105459, -1.601788], [-0.628434, 1.058759, -0.340964], [-0.082557, 0.188946, 1.491612], [-6.890503, -2.06899, -1.960068]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_Intra_Diels_alder_R_9_10', 'label': '5009017845242299296281_Intra_Diels_alder_R_9_10', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': ''}
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
