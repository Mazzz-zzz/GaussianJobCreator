import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0177'
logfile = 'conf/5009017845242299296281_0177.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863822, 0.7718203945763815, 1.163533622908849], [-2.2709622836291916, 0.7431123812655613, 1.1797556627389014], [-2.9970239643018997, 1.4656877041954783, 2.3681573973696826], [-2.270485399678259, 2.796025839123985, 2.7705593355752014], [-3.0948301649158365, 3.5679399740158084, 3.456931114902275], [-1.1948366276544546, 2.5317466231874195, 3.506790451340248], [-1.6927847351456848, 3.757652723316707, 1.2661660994459583], [-2.7310271787076, 3.7296051822020764, 0.29472543134830986], [-1.086022838358945, 4.944294898290884, 1.7447861166081318], [-0.5420771395216657, 2.7868089324098895, 0.8036977298898708], [-4.239914773111563, 1.7601618341371341, 1.997348852865152], [-3.023054582933851, 0.6702887991292501, 3.4298036615283656], [-2.6427825570536143, 1.3431159399560015, 0.04548369333676736], [-2.667789377892029, -0.5232706377056777, 1.1465664052870914], [-0.3710451618282828, 2.0617372872159105, 1.060259129110613], [-0.24552532002048744, 0.3056731502912612, 2.324090564665854], [1.5770424436171666, 0.0, 0.0], [2.2927181468939173, 1.391552724358058, 0.0], [2.341079856722324, 2.05985269289493, 1.4165023767064717], [1.145391641007038, 2.0331725558915963, 1.975951122830792], [3.202217992284646, 1.4422259961598325, 2.1984030628363267], [2.721555559149504, 3.3152059779287955, 1.2745358845394026], [1.6292323391939756, 2.2122558673105766, -0.8090479336198871], [3.5455868300943836, 1.2600392214310756, -0.4280914688619842], [1.9974224573334838, -0.6906780683055247, 1.0535722235492986], [1.9277183224308965, -0.6529932317206242, -1.1102241252095326], [-0.35014935725347407, -1.2838136616209443, 0.08241309473864848], [-0.4266843221927565, 0.4915335256355435, -1.1586058166012256], [-0.22149308850549337, 2.2648848865288254, 1.5525525511077085]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0177', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
