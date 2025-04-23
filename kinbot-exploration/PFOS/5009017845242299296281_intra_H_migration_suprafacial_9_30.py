import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_intra_H_migration_suprafacial_9_30'
logfile = '5009017845242299296281_intra_H_migration_suprafacial_9_30.log'

atom = ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H']
geom = [[0.740132, -2.881445, -5.478409], [1.431781, -3.378059, -4.163271], [2.204826, -2.210388, -3.482337], [1.355585, -1.612682, -2.279999], [2.018473, -0.322732, -1.787682], [3.305244, -0.250609, -2.069542], [1.499225, 0.782207, -2.333741], [1.297817, 0.126155, 0.255463], [1.957534, -0.175967, -0.46203], [2.202102, 1.060675, 1.255491], [0.297873, -0.036412, 0.589777], [0.095344, -1.395952, -2.667078], [1.255465, -2.523047, -1.305909], [3.386029, -2.640401, -3.033783], [2.527663, -1.250266, -4.351117], [2.263623, -4.382412, -4.423814], [0.534469, -3.929161, -3.323103], [-0.252669, -4.008921, -5.988928], [-0.749509, -3.639471, -7.421956], [-1.890256, -4.601892, -7.843253], [-1.578462, -5.87961, -7.691396], [-3.008686, -4.405191, -7.152761], [-2.219247, -4.45407, -9.121171], [0.254453, -3.684307, -8.300379], [-1.165535, -2.368318, -7.478507], [0.34433, -5.202608, -5.984473], [-1.27164, -4.146591, -5.139881], [0.093855, -1.74435, -5.292585], [1.65245, -2.600188, -6.406232], [0.760955, -0.2606, -0.494333]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_intra_H_migration_suprafacial_9_30', 'label': '5009017845242299296281_intra_H_migration_suprafacial_9_30', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': ''}
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
