import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0113'
logfile = 'conf/5009017845242299296281_0113.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863853, 0.621739478308218, -1.250182880316499], [-0.3976197158559584, -0.07566485901595763, -2.633910119820668], [-1.123388212146639, 0.4940238789133656, -3.902883790150536], [-1.2235971086413808, 2.0585875945462284, -3.8618506841213947], [-2.2148482487278947, 2.4243645240243836, -3.0680431770923215], [-0.08476421785157753, 2.5840652181930435, -3.4193205106248867], [-1.546098718177794, 2.781561421086142, -5.563310318993654], [-2.0148297746218793, 4.114100839171646, -5.398808331822813], [-0.44471162371596595, 2.4228985185350953, -6.37796617402361], [-2.760362629999281, 1.8660327814347824, -5.973082319300063], [-0.43435958530525626, 0.14440355649341152, -4.98547176997021], [-2.3533195674001606, 0.0018245094686658277, -3.976629593905572], [0.9186702629217383, 0.059262854653436996, -2.8183083047232307], [-0.6863393184328432, -1.3666682794365153, -2.5232286857846646], [-0.25604457595343316, 1.8735740976390431, -1.3916701657561947], [-2.00760247718745, 0.6427130616946718, -1.0543092166280599], [1.5770424436171635, 0.0, 0.0], [2.2927181468939093, 1.391552724358059, 0.0], [2.341079856722318, 2.0598526928949337, 1.4165023767064735], [1.1453916410070293, 2.0331725558915896, 1.9759511228307978], [3.2022179922846408, 1.4422259961598314, 2.1984030628363267], [2.7215555591494933, 3.3152059779288026, 1.2745358845394108], [1.6292323391939667, 2.212255867310583, -0.8090479336198841], [3.545586830094378, 1.260039221431088, -0.4280914688619807], [1.997422457333488, -0.6906780683055223, 1.0535722235492966], [1.9277183224308938, -0.6529932317206203, -1.110224125209535], [-0.3501493572534765, 0.7132786644586344, 1.070608697319935], [-0.42668432219275737, -1.2491488329668539, 0.15362238828850275], [-2.7315131263772168, 1.0312081612344972, -5.485027049536264]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0113', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
