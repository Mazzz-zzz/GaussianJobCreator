import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_well_mp2'
logfile = '5009017845242299296281_well_mp2.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[-1.111068, 0.248152, -0.510411], [-0.388429, -0.421649, 0.697867], [1.066889, -0.967414, 0.427785], [2.139386, 0.053534, -0.090862], [3.599565, -0.439822, 0.199073], [3.910989, -0.203548, 1.461488], [3.705578, -1.742562, -0.046955], [4.85896, 0.435209, -0.88255], [4.801847, -0.136662, -2.18337], [4.73504, 1.821057, -0.619062], [6.166234, -0.071514, -0.165048], [1.996422, 0.196595, -1.405429], [1.971419, 1.225542, 0.507982], [1.486068, -1.441224, 1.604484], [0.984472, -1.972254, -0.435814], [-1.096353, -1.489387, 1.068616], [-0.341761, 0.438546, 1.709281], [-2.553213, 0.823556, -0.234384], [-3.622558, -0.167592, 0.333716], [-5.094565, 0.326725, 0.123267], [-5.215766, 1.58217, 0.513822], [-5.44731, 0.227592, -1.141697], [-5.89766, -0.428595, 0.848083], [-3.426963, -0.301817, 1.64217], [-3.505561, -1.353272, -0.258494], [-2.443432, 1.850106, 0.601028], [-2.982985, 1.273988, -1.415222], [-0.385555, 1.295394, -0.903462], [-1.184449, -0.626571, -1.508162], [6.519252, -0.849088, -0.619712]])

kwargs = {'method': 'mp2', 'basis': '6-31G', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_well_mp2', 'label': '5009017845242299296281_well_mp2', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
Gaussian.command = 'g16 < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy()  # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    freq = reader_gauss.read_freq(logfile, ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'])
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
            freq = reader_gauss.read_freq(logfile, ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'])
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
