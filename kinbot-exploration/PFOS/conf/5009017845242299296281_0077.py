import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0077'
logfile = 'conf/5009017845242299296281_0077.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863819, -1.3935598728845984, 0.08664925740765084], [-0.39761971585595535, -2.2432006455416347, 1.3824827499919197], [-0.7491833517666578, -1.5977521800128374, 2.768620777958024], [0.37513150925004024, -0.6189636534602634, 3.255934098715886], [0.8222268376054811, 0.09417028242405033, 2.237087839494269], [-0.09608992658046997, 0.19684114923706436, 4.194716780097906], [1.8294372022557688, -1.5425873314445653, 3.9999620443421735], [2.942173082714334, -0.6581403464427761, 4.0447213882977495], [1.3352845675176521, -2.2532055988734037, 5.120761590749029], [2.0568917557993647, -2.5820463860490737, 2.838785085527545], [-1.882037352830636, -0.911344475583573, 2.6482969312796576], [-0.8963077255402113, -2.547907795232932, 3.682944859997027], [-1.130548978421929, -3.3520534153271306, 1.247992163869569], [0.8878652548597811, -2.5746029018458154, 1.382089214148236], [-0.256044575953427, -2.1420087660532867, -0.9267276815498231], [-2.0076024771874446, -1.234415095891311, -0.02945123045763413], [1.577042443617165, 0.0, 0.0], [2.2927181468939164, 1.3915527243580552, 0.0], [1.6005215470082557, 2.440721904563901, -0.9357086002340301], [1.3760692761371665, 1.9145350254105769, -2.125714308102146], [0.4611709133627573, 2.851946778983254, -0.41882897761344917], [2.4048649076935025, 3.478080011182911, -1.068930340330689], [3.5400592233304637, 1.2295174124846067, -0.43152105569274785], [2.3088468039522523, 1.896094738758371, 1.2310220414904698], [1.9974224573334847, -0.6906780683055274, 1.0535722235492977], [1.9277183224308956, -0.6529932317206254, -1.1102241252095333], [-0.3501493572534722, 0.5705349971623088, -1.1530217920585821], [-0.42668432219275876, 0.757615307331304, 1.0049834283127257], [1.239325834061783, -2.7103870890369044, 2.3377092287294263]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0077', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
