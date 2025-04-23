import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0497'
logfile = 'conf/5009017845242299296281_0497.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863834, -1.3935598728845962, 0.08664925740765816], [-2.270962283629192, -1.3932545648232157, 0.05367636867328109], [-2.997023964301902, -2.783728318379923, 0.08524408716294642], [-4.4759542789311215, -2.655862593215632, 0.591209145048298], [-5.168854707529478, -3.720532137165071, 0.2266961805165777], [-4.503853254513124, -2.5414503697214306, 1.9159851754887638], [-5.326480128927349, -1.1401075478729004, -0.1164318318045449], [-4.9901954846991545, -1.0442971716297287, -1.4948962011918954], [-6.65347461560192, -1.1448527563795778, 0.37795656508552794], [-4.5262555985949, -0.039328806966301096, 0.6763288473611787], [-3.015076885801345, -3.287350694771378, -1.1456215214756071], [-2.3510956819310227, -3.610181153269013, 0.8977600360686799], [-2.6427825570536148, -0.7109480038655783, 1.140430677561336], [-2.6677893778920336, -0.7313203152515878, -1.0264488679511374], [-0.3710451618282827, -1.949079984012109, 1.2553873021032933], [-0.2455253200204909, -2.16555804484198, -0.8973245689258723], [1.5770424436171664, 0.0, 0.0], [2.29271814689392, 1.3915527243580559, 0.0], [3.782335574419715, 1.318614735245459, -0.4807937764724458], [4.419592497958774, 0.3492507324843921, 0.1497631852713449], [3.8473329977166344, 1.1104432852891963, -1.779574085222885], [4.368010411051583, 2.4684340853352005, -0.20560554420872534], [2.293896732220221, 1.8704189044736026, 1.240568989312628], [1.6494649440008824, 2.2352148943143333, -0.8029305726284919], [1.9974224573334851, -0.6906780683055256, 1.0535722235493008], [1.9277183224308918, -0.6529932317206285, -1.1102241252095326], [-0.3501493572534752, 0.5705349971623087, -1.1530217920585795], [-0.4266843221927528, 0.7576153073313044, 1.0049834283127241], [-4.169402957650365, -0.4084543011389569, 1.4963074744599043]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0497', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
