import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0352'
logfile = 'conf/5009017845242299296281_0352.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863831, 0.7718203945763825, 1.1635336229088504], [-0.3976197158559554, 2.3188655045575883, 1.2514273698287528], [-0.7491833517666567, 3.1965720171635055, -0.0006164121359225181], [-0.9314389667275255, 4.706813247194561, 0.3808799534188789], [-0.008749002261886999, 5.063027021493548, 1.2573826171633007], [-0.8331735067029749, 5.473001427583927, -0.7018017228603524], [-2.614049859498326, 5.035074094275036, 1.1442466614419928], [-3.5783758997940347, 5.064491441695677, 0.09943045009746962], [-2.712669687109959, 4.211324697515599, 2.291882419245102], [-2.357101367856443, 6.509078486159019, 1.6360099735736495], [0.2407930719723934, 3.1018030114350617, -0.8837681932585539], [-1.8795656402123513, 2.7679819246717727, -0.5474035883527788], [-1.1305489784219307, 2.75681962529852, 2.278967330580905], [0.8878652548597803, 2.4842258206717465, 1.5386269105814976], [-0.2560445759534309, 0.2684346684142367, 2.3183978473060143], [-2.007602477187446, 0.5917020341966337, 1.0837604470856974], [1.5770424436171653, 0.0, 0.0], [2.2927181468939146, 1.3915527243580554, 0.0], [3.782335574419715, 1.3186147352454685, -0.48079377647244503], [4.4195924979587735, 0.34925073248439353, 0.14976318527135224], [3.847332997716637, 1.1104432852892046, -1.779574085222877], [4.368010411051579, 2.4684340853352094, -0.20560554420871324], [2.293896732220217, 1.8704189044736022, 1.2405689893126362], [1.6494649440008788, 2.23521489431434, -0.8029305726284843], [1.9974224573334813, -0.6906780683055324, 1.0535722235492966], [1.9277183224308936, -0.6529932317206242, -1.1102241252095362], [-0.3501493572534752, -1.283813661620945, 0.08241309473864737], [-0.4266843221927584, 0.49153352563554864, -1.1586058166012276], [-2.8140034156581577, 6.668305929388894, 2.47377457460014]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0352', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
