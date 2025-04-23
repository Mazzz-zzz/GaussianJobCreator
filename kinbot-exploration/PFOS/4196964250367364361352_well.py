import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '4196964250367364361352_well'
logfile = '4196964250367364361352_well.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F'], positions=[[0.029811, 0.077906, -0.035169], [-0.691473, -1.294731, 0.158236], [-2.268621, -1.259468, 0.148804], [-2.986577, -0.728806, -1.144108], [-4.448908, -1.281274, -1.263868], [-4.425072, -2.527998, -1.722724], [-5.049465, -1.253671, -0.08471], [-3.055232, 0.597941, -1.092506], [-2.305626, -1.101698, -2.220917], [-2.643406, -2.529521, 0.330264], [-2.673812, -0.536116, 1.185604], [-0.363371, -1.782965, 1.352356], [-0.272763, -2.136042, -0.784041], [1.589669, -0.06159, -0.178601], [2.328409, 1.295276, 0.090139], [3.79174, 1.321567, -0.464947], [4.444971, 0.242144, -0.073635], [3.792388, 1.377048, -1.780944], [4.404736, 2.389999, 0.006482], [2.385089, 1.500772, 1.401244], [1.654434, 2.290294, -0.483316], [2.045481, -0.965254, 0.683495], [1.863299, -0.468264, -1.416204], [-0.431907, 0.650374, -1.14422], [-0.248448, 0.853877, 1.009781]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'chk': '4196964250367364361352_well', 'label': '4196964250367364361352_well', 'Symm': 'None', 'mult': 2, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
Gaussian.command = 'g16 < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy()  # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    freq = reader_gauss.read_freq(logfile, ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F'])
    zpe = reader_gauss.read_zpe(logfile)
    db.write(mol, name=label, data={'energy': e, 'frequencies': np.asarray(freq),
                                     'zpe': zpe, 'status': 'normal'})

except RuntimeError:
    for i in range(3):
        try:
            iowait(logfile, 'gauss')
            mol.positions = reader_gauss.read_geom(logfile, mol)
            kwargs = reader_gauss.correct_kwargs(logfile, kwargs)
            mol.calc = Gaussian(**kwargs)
            e = mol.get_potential_energy()  # use the Gaussian optimizer
            iowait(logfile, 'gauss')
            mol.positions = reader_gauss.read_geom(logfile, mol)
            freq = reader_gauss.read_freq(logfile, ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F'])
            zpe = reader_gauss.read_zpe(logfile)
            db.write(mol, name=label, data={'energy': e,
                                             'frequencies': np.asarray(freq),
                                             'zpe': zpe, 'status': 'normal'})
        except RuntimeError:
            if i == 2:
                db.write(mol, name=label, data={'status': 'error'})
            pass
        else:
            break

with open(logfile, 'a') as f:
    f.write('done\n')
