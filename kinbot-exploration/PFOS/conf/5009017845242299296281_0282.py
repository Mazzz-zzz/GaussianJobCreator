import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0282'
logfile = 'conf/5009017845242299296281_0282.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863829, 0.6217394783082163, -1.2501828803165014], [-0.3466020415139031, 2.127818130564327, -1.5663863369811213], [1.1624535490467067, 2.4914860914100183, -1.7945884035232582], [1.9243606525366, 2.688842566201291, -0.43807332933071275], [3.2268445905155723, 2.58295218452213, -0.6351133289542098], [1.6500667572714027, 3.883782491250017, 0.07746665466450789], [1.4242023469536254, 1.404868453778162, 0.8358636703276711], [2.4092408723075187, 1.3886830619286343, 1.8614699029868225], [0.03808133731720843, 1.5760028148920204, 1.069825950917271], [1.6020493994449863, 0.11866616569049313, -0.05553788332095845], [1.7462981079693125, 1.502568076923133, -2.4655267346281247], [1.2535978601361857, 3.618508316135838, -2.488809853075985], [-1.0117971004258968, 2.4109962356139936, -2.689842246735281], [-0.8215282967123797, 2.8830315798638466, -0.5833314461381539], [-2.0119541879597174, 0.6042709716797826, -1.0485616399675899], [-0.40994706586377827, -0.11587296658230002, -2.3181096973944286], [1.5770424436171664, 0.0, 0.0], [2.2927181468939173, 1.3915527243580568, 0.0], [1.60052154700825, 2.4407219045639, -0.9357086002340265], [1.3760692761371622, 1.9145350254105744, -2.1257143081021432], [0.4611709133627542, 2.8519467789832547, -0.41882897761344384], [2.404864907693497, 3.4780800111829104, -1.0689303403306853], [3.5400592233304637, 1.2295174124846042, -0.43152105569274857], [2.3088468039522545, 1.8960947387583746, 1.2310220414904698], [1.997422457333487, -0.6906780683055267, 1.0535722235492981], [1.9277183224308947, -0.6529932317206247, -1.1102241252095362], [-0.35014935725347346, 0.713278664458634, 1.0706086973199342], [-0.42668432219275687, -1.2491488329668508, 0.1536223882885001], [0.926447763097055, -0.538193910108853, 0.16368724407294347]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0282', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
