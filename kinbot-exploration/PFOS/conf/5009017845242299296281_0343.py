import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0343'
logfile = 'conf/5009017845242299296281_0343.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863827, 0.7718203945763841, 1.1635336229088504], [-0.3976197158559555, 2.318865504557592, 1.251427369828748], [-0.7491833517666581, 3.1965720171635046, -0.0006164121359296159], [-0.9314389667275264, 4.706813247194561, 0.3808799534188684], [-0.008749002261888154, 5.06302702149355, 1.2573826171632891], [-0.8331735067029783, 5.473001427583925, -0.7018017228603642], [-2.6140498594983286, 5.035074094275036, 1.1442466614419813], [-2.8987158918487608, 3.986054214870736, 2.061281857304407], [-2.653689516643581, 6.411227428777662, 1.4758902489975207], [-3.4932515493432246, 4.830132937376362, -0.14631373644411508], [0.240793071972392, 3.1018030114350594, -0.8837681932585609], [-1.8795656402123533, 2.7679819246717687, -0.5474035883527844], [-1.1305489784219305, 2.756819625298526, 2.278967330580899], [0.8878652548597799, 2.4842258206717514, 1.5386269105814927], [-0.25604457595342867, 0.2684346684142433, 2.3183978473060174], [-2.007602477187445, 0.5917020341966351, 1.0837604470856999], [1.5770424436171653, 0.0, 0.0], [2.292718146893916, 1.3915527243580565, 0.0], [2.3410798567223283, 2.0598526928949337, 1.4165023767064735], [1.1453916410070368, 2.033172555891595, 1.9759511228307916], [3.2022179922846465, 1.4422259961598296, 2.1984030628363254], [2.721555559149504, 3.3152059779288, 1.2745358845394028], [1.6292323391939747, 2.212255867310577, -0.8090479336198863], [3.545586830094382, 1.2600392214310736, -0.4280914688619836], [1.9974224573334811, -0.69067806830553, 1.0535722235492984], [1.9277183224308951, -0.6529932317206266, -1.110224125209531], [-0.35014935725347446, -1.2838136616209457, 0.08241309473865084], [-0.42668432219275737, 0.4915335256355448, -1.1586058166012259], [-2.9652438374406374, 4.988680915744008, -0.9413186304057337]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0343', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
