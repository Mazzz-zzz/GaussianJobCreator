import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0377'
logfile = 'conf/5009017845242299296281_0377.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863807, 0.7718203945763812, 1.163533622908853], [-0.3976197158559554, 2.318865504557588, 1.251427369828754], [-0.7491833517666567, 3.1965720171635055, -0.0006164121359210986], [0.37513150925003946, 3.129203469266082, -1.0919288014421193], [1.3770046036719092, 3.9220036113121046, -0.7539640391926383], [0.8207913513067462, 1.882785072952196, -1.2212676086939316], [-0.25282578633505465, 3.6732379270234503, -2.7745931990000794], [0.8709723088773929, 3.9682319665744465, -3.5946464689992017], [-1.2861204129680759, 2.7755913864028683, -3.1377593790488643], [-0.9289189133715272, 5.027138424149858, -2.338251409500251], [-1.882037352830636, 2.7491646570443433, -0.5349009981858397], [-0.8963077255402094, 4.463477707111216, 0.36508044717361227], [-1.1305489784219294, 2.7568196252985206, 2.278967330580905], [0.8878652548597796, 2.484225820671747, 1.5386269105814978], [-0.25604457595342545, 0.26843466841423697, 2.318397847306017], [-2.007602477187444, 0.5917020341966325, 1.0837604470856996], [1.5770424436171635, 0.0, 0.0], [2.2927181468939155, 1.3915527243580534, 0.0], [2.3410798567223288, 2.05985269289493, 1.4165023767064695], [1.1453916410070415, 2.0331725558915985, 1.9759511228307924], [3.202217992284651, 1.4422259961598263, 2.1984030628363227], [2.721555559149507, 3.315205977928798, 1.2745358845393986], [1.6292323391939734, 2.2122558673105774, -0.8090479336198887], [3.5455868300943805, 1.2600392214310732, -0.428091468861989], [1.997422457333483, -0.6906780683055296, 1.0535722235492968], [1.9277183224308931, -0.6529932317206221, -1.1102241252095353], [-0.3501493572534767, -1.283813661620945, 0.08241309473864623], [-0.42668432219275854, 0.49153352563554764, -1.1586058166012254], [-1.1676786210542363, 4.99684909362526, -1.4012141929823398]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0377', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
