import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0240'
logfile = 'conf/5009017845242299296281_0240.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863867, -1.3935598728845966, 0.08664925740765446], [-0.34660204151390483, -2.4204394252486705, -1.0595513872112494], [-0.987392984244527, -3.8498102687678504, -0.9708648395635502], [-1.0274356123761075, -4.556709010365346, -2.370248959832962], [-1.1749864386042315, -5.860276414715168, -2.210741219979811], [-2.0336420090655967, -4.083012157772831, -3.0996734039381115], [0.5505558709515881, -4.26641364969085, -3.343560222087291], [0.5065706167013714, -2.9499985957706922, -3.8797949757455483], [1.6165458359829274, -4.767649857467604, -2.5575675821000754], [0.2920999651411268, -5.29495502384334, -4.50796814726524], [-0.26149276674437044, -4.595762031236107, -0.14284527771145913], [-2.230877915023902, -3.7580110466497003, -0.5172703827718951], [-0.7784534903451407, -1.8458422644301347, -2.1856818928204538], [0.973158224737943, -2.549721586648366, -1.1213428320614764], [-2.0119541879597196, -1.210216503485696, 0.0009668077396009834], [-0.4099470658637832, -1.949605403411477, 1.259403781369362], [1.5770424436171624, 0.0, 0.0], [2.2927181468939155, 1.3915527243580552, 0.0], [2.3410798567223208, 2.059852692894929, 1.4165023767064744], [1.1453916410070348, 2.033172555891597, 1.9759511228307898], [3.2022179922846434, 1.4422259961598223, 2.198403062836329], [2.7215555591495066, 3.3152059779287915, 1.2745358845394072], [1.6292323391939858, 2.212255867310575, -0.8090479336198897], [3.5455868300943827, 1.2600392214310676, -0.4280914688619775], [1.9974224573334805, -0.6906780683055251, 1.053572223549304], [1.9277183224308958, -0.6529932317206331, -1.110224125209526], [-0.3501493572534783, 0.5705349971623067, -1.1530217920585857], [-0.42668432219276287, 0.7576153073313088, 1.0049834283127244], [1.123239438852584, -5.706186785466734, -4.783754296603381]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0240', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
