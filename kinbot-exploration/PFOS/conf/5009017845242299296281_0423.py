import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0423'
logfile = 'conf/5009017845242299296281_0423.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863838, 0.7718203945763857, 1.1635336229088462], [-0.39761971585595596, 2.3188655045575928, 1.2514273698287444], [1.087941409756326, 2.759867566386747, 1.496793764948379], [1.3412747433710004, 4.235018190162302, 1.0284583955214368], [0.3135944378492583, 4.99529163351123, 1.3636589683999776], [2.4483457429192925, 4.714817055771962, 1.5880000749983478], [1.5638064896882018, 4.347671144906446, -0.8315018902866818], [1.387731031182455, 5.704567054301786, -1.2191382298991764], [2.702177123588236, 3.5708891462475, -1.157199292919017], [0.28901432478710026, 3.5296396286834955, -1.263094295898448], [1.3507067731889253, 2.680176727710471, 2.798175832796711], [1.9101736356828996, 1.9646745303205164, 0.8244605116218734], [-0.7843720924895173, 2.815878982142089, 0.07319662459687076], [-1.1581983763610544, 2.837834268725447, 2.207717463520171], [-0.2560445759534288, 0.26843466841424485, 2.318397847306013], [-2.007602477187449, 0.591702034196638, 1.0837604470856914], [1.5770424436171646, 0.0, 0.0], [2.292718146893915, 1.391552724358055, 0.0], [1.6005215470082501, 2.4407219045638984, -0.9357086002340297], [1.376069276137164, 1.914535025410573, -2.125714308102146], [0.46117091336275506, 2.851946778983252, -0.41882897761345506], [2.4048649076934967, 3.4780800111829113, -1.068930340330685], [3.540059223330463, 1.229517412484606, -0.43152105569274535], [2.308846803952248, 1.896094738758374, 1.23102204149047], [1.9974224573334827, -0.6906780683055262, 1.0535722235492995], [1.9277183224308927, -0.6529932317206277, -1.1102241252095313], [-0.3501493572534762, -1.2838136616209452, 0.08241309473865081], [-0.42668432219275937, 0.49153352563554376, -1.1586058166012299], [0.021434877808575684, 2.926374443114548, -0.5556801723957191]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0423', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
