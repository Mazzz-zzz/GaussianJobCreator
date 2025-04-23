import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0060'
logfile = 'conf/5009017845242299296281_0060.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863802, -1.3935598728845995, 0.0866492574076509], [-0.3976197158559543, -2.243200645541638, 1.382482749991916], [-0.7491833517666566, -1.5977521800128405, 2.768620777958024], [0.37513150925004174, -0.6189636534602692, 3.2559340987158847], [1.3770046036719104, -1.3080497941752989, 3.7735367807269125], [0.8207913513067474, 0.11625623747192053, 2.2411735073897048], [-0.25282578633505265, 0.5662492319898811, 4.568413958446837], [-1.0611560354132719, -0.16726554854875333, 5.480000094237909], [0.8456628214350043, 1.3726180086267916, 4.953704688042676], [-1.1898482641172174, 1.4421209781482605, 3.654422607163112], [-1.8820373528306358, -0.9113444755835776, 2.648296931279656], [-0.8963077255402084, -2.5479077952329363, 3.6829448599970265], [-1.1305489784219294, -3.3520534153271324, 1.2479921638695661], [0.8878652548597796, -2.5746029018458176, 1.382089214148231], [-0.2560445759534265, -2.1420087660532863, -0.9267276815498242], [-2.0076024771874446, -1.2344150958913134, -0.02945123045763419], [1.5770424436171646, 0.0, 0.0], [2.2927181468939146, 1.3915527243580577, 0.0], [1.6005215470082486, 2.4407219045638993, -0.9357086002340318], [1.3760692761371645, 1.914535025410578, -2.125714308102144], [0.46117091336275085, 2.851946778983253, -0.4188289776134476], [2.4048649076934963, 3.478080011182911, -1.06893034033069], [3.540059223330461, 1.2295174124846096, -0.43152105569274846], [2.308846803952253, 1.8960947387583746, 1.231022041490469], [1.9974224573334838, -0.6906780683055267, 1.053572223549296], [1.927718322430895, -0.652993231720622, -1.1102241252095333], [-0.3501493572534758, 0.5705349971623108, -1.15302179205858], [-0.4266843221927588, 0.7576153073313024, 1.0049834283127272], [-0.8961744767104263, 1.3996910745483768, 2.7335980871992622]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0060', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
