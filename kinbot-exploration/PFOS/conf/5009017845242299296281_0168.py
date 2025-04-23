import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0168'
logfile = 'conf/5009017845242299296281_0168.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863838, 0.7718203945763857, 1.1635336229088462], [-0.39761971585595596, 2.3188655045575928, 1.2514273698287444], [1.087941409756326, 2.759867566386747, 1.496793764948379], [1.3412747433710004, 4.235018190162302, 1.0284583955214368], [0.3135944378492583, 4.99529163351123, 1.3636589683999776], [2.4483457429192925, 4.714817055771962, 1.5880000749983478], [1.5638064896882018, 4.347671144906446, -0.8315018902866818], [1.387731031182455, 5.704567054301786, -1.2191382298991764], [2.702177123588236, 3.5708891462475, -1.157199292919017], [0.28901432478710026, 3.5296396286834955, -1.263094295898448], [1.3507067731889253, 2.680176727710471, 2.798175832796711], [1.9101736356828996, 1.9646745303205164, 0.8244605116218734], [-0.7843720924895173, 2.815878982142089, 0.07319662459687076], [-1.1581983763610544, 2.837834268725447, 2.207717463520171], [-0.2560445759534288, 0.26843466841424485, 2.318397847306013], [-2.007602477187449, 0.591702034196638, 1.0837604470856914], [1.5770424436171646, 0.0, 0.0], [2.292718146893915, 1.391552724358055, 0.0], [2.3410798567223217, 2.0598526928949283, 1.4165023767064744], [1.1453916410070333, 2.0331725558915923, 1.9759511228307929], [3.20221799228464, 1.4422259961598298, 2.1984030628363285], [2.721555559149502, 3.3152059779287955, 1.2745358845394055], [1.6292323391939756, 2.212255867310575, -0.8090479336198863], [3.5455868300943805, 1.2600392214310734, -0.4280914688619809], [1.997422457333482, -0.690678068305526, 1.0535722235492997], [1.9277183224308925, -0.6529932317206276, -1.1102241252095313], [-0.3501493572534762, -1.2838136616209452, 0.08241309473865081], [-0.42668432219275937, 0.49153352563554376, -1.1586058166012299], [-0.4465738568915164, 4.128826514682334, -1.452428628014321]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0168', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
