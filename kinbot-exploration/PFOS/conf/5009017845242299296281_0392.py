import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0392'
logfile = 'conf/5009017845242299296281_0392.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863874, 0.621739478308216, -1.2501828803164978], [-0.3466020415139043, 2.127818130564331, -1.5663863369811153], [-0.9873929842445258, 2.7656987490870812, -2.848601072721362], [-2.462172495902625, 2.278219176447387, -3.065336434982642], [-2.461597984833636, 1.071414784733829, -3.603968752203674], [-3.1141435203244816, 2.2457738046910265, -1.9065521057311172], [-3.4141816509334397, 3.422830080870185, -4.207793338908972], [-3.7843994872934696, 4.580728877090411, -3.469795954657076], [-2.7102711120366747, 3.4564154765019377, -5.436100607501947], [-4.692190264502985, 2.5260059550865663, -4.4150205110257374], [-0.9969507321948232, 4.088639050760538, -2.711746550807541], [-0.2788854647228134, 2.4314306036058904, -3.9194473283466773], [-0.7784534903451381, 2.815777175989243, -0.5057053459652454], [0.9731582247379434, 2.24597217224102, -1.6474522505842997], [-2.011954187959716, 0.6042709716797846, -1.048561639967588], [-0.409947065863779, -0.1158729665822953, -2.318109697394427], [1.5770424436171655, 0.0, 0.0], [2.292718146893918, 1.3915527243580603, 0.0], [2.3410798567223265, 2.05985269289493, 1.4165023767064764], [1.1453916410070384, 2.033172555891593, 1.9759511228307916], [3.2022179922846483, 1.4422259961598307, 2.1984030628363223], [2.7215555591495004, 3.315205977928799, 1.2745358845394024], [1.629232339193972, 2.2122558673105788, -0.8090479336198849], [3.545586830094379, 1.2600392214310807, -0.42809146886198507], [1.9974224573334882, -0.6906780683055254, 1.0535722235492977], [1.9277183224308918, -0.652993231720621, -1.11022412520953], [-0.3501493572534726, 0.7132786644586342, 1.0706086973199367], [-0.4266843221927558, -1.2491488329668483, 0.15362238828850205], [-5.38940782022126, 2.78658795577981, -3.797001832012762]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0392', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
