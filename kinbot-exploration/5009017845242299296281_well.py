import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '5009017845242299296281_well'
logfile = '5009017845242299296281_well.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[-1.046297, 0.131055, -0.504729], [-0.37462, -0.322492, 0.771945], [1.01986, -0.91147, 0.606663], [2.054013, 0.005959, -0.007744], [3.441131, -0.582763, 0.151439], [3.771559, -0.71443, 1.468577], [3.491916, -1.841356, -0.374827], [4.607156, 0.521506, -0.716157], [4.478503, 0.397087, -2.150039], [4.681339, 1.81903, -0.077732], [5.996988, -0.243994, -0.373882], [1.776261, 0.180928, -1.327455], [1.976059, 1.231123, 0.569154], [1.442032, -1.30086, 1.847147], [0.948385, -2.066437, -0.11998], [-1.142286, -1.285461, 1.36393], [-0.340204, 0.709405, 1.65896], [-2.441772, 0.720157, -0.346165], [-3.472246, -0.203554, 0.259891], [-4.861591, 0.390794, 0.101631], [-4.966157, 1.603318, 0.68605], [-5.217584, 0.534199, -1.19718], [-5.798549, -0.405263, 0.673672], [-3.211083, -0.40443, 1.577961], [-3.410395, -1.430283, -0.331494], [-2.372135, 1.868193, 0.385673], [-2.857798, 1.112886, -1.587309], [-0.276338, 1.090518, -1.093632], [-1.080025, -0.904988, -1.392241], [6.428372, -0.513181, -1.204099]])

kwargs = {'method': 'bmk', 'basis': '6-31++G(2df,p)', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_well', 'label': '5009017845242299296281_well', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
