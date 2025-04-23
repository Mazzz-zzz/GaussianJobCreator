import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0048'
logfile = 'conf/5009017845242299296281_0048.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586383, -1.393559872884599, 0.08664925740765211], [-0.39761971585595723, -2.2432006455416356, 1.3824827499919201], [-0.7491833517666586, -1.5977521800128376, 2.7686207779580245], [0.37513150925003946, -0.6189636534602634, 3.255934098715886], [0.8222268376054808, 0.09417028242405236, 2.2370878394942695], [-0.09608992658047184, 0.19684114923706436, 4.194716780097906], [1.8294372022557648, -1.542587331444566, 3.9999620443421753], [1.4673342643132088, -1.961144854839958, 5.309970467884738], [2.3207428086512274, -2.412711909526456, 2.996538508944387], [2.8319974601065314, -0.33339637242824155, 4.115169612758457], [-1.8820373528306378, -0.911344475583573, 2.6482969312796576], [-0.8963077255402113, -2.5479077952329288, 3.6829448599970296], [-1.1305489784219303, -3.352053415327129, 1.2479921638695717], [0.8878652548597781, -2.5746029018458167, 1.3820892141482366], [-0.2560445759534286, -2.142008766053288, -0.9267276815498215], [-2.007602477187445, -1.2344150958913112, -0.029451230457635234], [1.5770424436171648, 0.0, 0.0], [2.2927181468939164, 1.3915527243580552, 0.0], [3.7823355744197156, 1.318614735245459, -0.4807937764724444], [4.419592497958777, 0.3492507324843942, 0.14976318527135551], [3.847332997716638, 1.1104432852891963, -1.7795740852228772], [4.368010411051583, 2.4684340853352054, -0.20560554420871502], [2.2938967322202197, 1.8704189044736035, 1.24056898931263], [1.6494649440008817, 2.2352148943143373, -0.8029305726284852], [1.9974224573334833, -0.6906780683055286, 1.0535722235492977], [1.9277183224308956, -0.6529932317206258, -1.1102241252095326], [-0.35014935725347224, 0.5705349971623084, -1.1530217920585826], [-0.42668432219275865, 0.7576153073313039, 1.0049834283127255], [2.6044076003362244, 0.354134980050604, 3.4737135489595845]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0048', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
