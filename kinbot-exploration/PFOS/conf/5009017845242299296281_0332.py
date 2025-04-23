import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0332'
logfile = 'conf/5009017845242299296281_0332.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863802, -1.3935598728845995, 0.08664925740765339], [-2.2709622836291903, -1.3932545648232186, 0.053676368673281204], [-2.997023964301899, -2.783728318379926, 0.08524408716294156], [-4.475954278931118, -2.655862593215639, 0.5912091450482946], [-4.492909775028066, -2.559896350243985, 1.909163975245072], [-5.0554016941443125, -1.5823291229684187, 0.06142935889679468], [-5.510548804511903, -4.144001785101999, 0.10467112564480845], [-5.849703873439191, -4.019393035774588, -1.2707842559005142], [-4.883527770796265, -5.282330859557806, 0.6671729595749599], [-6.787288074297166, -3.8263097931877783, 0.9704374134331277], [-3.015076885801341, -3.2873506947713813, -1.1456215214756134], [-2.3510956819310165, -3.6101811532690196, 0.8977600360686747], [-2.642782557053612, -0.7109480038655858, 1.1404306775613369], [-2.6677893778920327, -0.7313203152515898, -1.0264488679511357], [-0.371045161828278, -1.9490799840121131, 1.2553873021032886], [-0.2455253200204906, -2.16555804484198, -0.897324568925878], [1.5770424436171668, 0.0, 0.0], [2.292718146893917, 1.391552724358057, 0.0], [1.6005215470082492, 2.4407219045639015, -0.9357086002340242], [1.376069276137161, 1.9145350254105797, -2.12571430810214], [0.46117091336275595, 2.851946778983252, -0.41882897761344295], [2.404864907693498, 3.47808001118291, -1.0689303403306811], [3.540059223330462, 1.229517412484609, -0.43152105569274973], [2.3088468039522514, 1.896094738758374, 1.2310220414904705], [1.9974224573334867, -0.6906780683055253, 1.0535722235492988], [1.9277183224308942, -0.6529932317206266, -1.1102241252095368], [-0.35014935725347723, 0.57053499716231, -1.153021792058576], [-0.4266843221927535, 0.7576153073313006, 1.0049834283127284], [-7.178525391566183, -4.645289181981107, 1.3053530879003867]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0332', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
