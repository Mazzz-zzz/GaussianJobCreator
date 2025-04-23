import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0325'
logfile = 'conf/5009017845242299296281_0325.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863834, -1.3935598728845966, 0.08664925740765321], [-2.2709622836291934, -1.393254564823214, 0.05367636867327359], [-2.9970239643019023, -2.783728318379921, 0.08524408716293647], [-2.270485399678261, -3.79738768686225, 1.036149738531478], [-3.0948301649158383, -4.777760151646133, 1.3614610992245497], [-1.1948366276544569, -4.302842928203059, 0.4391616659556528], [-1.692784735145687, -2.975358369189201, 2.6211396672690723], [-1.4188572202552054, -3.9932116420477124, 3.575796380086378], [-0.7700728935849265, -1.9676476318325513, 2.2490124055592187], [-3.042335784893112, -2.2646639579155283, 3.0137820589169775], [-4.239914773111565, -2.609835763869493, 0.5256704367019975], [-3.023054582933852, -3.305441500441073, -1.134414702846084], [-2.6427825570536165, -0.710948003865578, 1.1404306775613289], [-2.6677893778920323, -0.7313203152515865, -1.0264488679511443], [-0.37104516182828423, -1.9490799840121085, 1.2553873021032882], [-0.245525320020489, -2.1655580448419807, -0.8973245689258739], [1.577042443617166, 0.0, 0.0], [2.292718146893919, 1.3915527243580554, 0.0], [2.3410798567223257, 2.0598526928949306, 1.416502376706471], [1.1453916410070368, 2.0331725558916003, 1.9759511228307898], [3.2022179922846465, 1.4422259961598263, 2.1984030628363267], [2.721555559149508, 3.3152059779287963, 1.2745358845394028], [1.6292323391939816, 2.212255867310575, -0.8090479336198881], [3.5455868300943854, 1.260039221431067, -0.428091468861979], [1.9974224573334807, -0.6906780683055227, 1.0535722235493017], [1.9277183224308954, -0.6529932317206283, -1.1102241252095324], [-0.3501493572534726, 0.5705349971623086, -1.1530217920585843], [-0.42668432219275587, 0.7576153073313057, 1.0049834283127221], [-2.8634284385792075, -1.4077394958496028, 3.425649773655509]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0325', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
