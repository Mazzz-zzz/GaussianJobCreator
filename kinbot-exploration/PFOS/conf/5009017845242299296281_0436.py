import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0436'
logfile = 'conf/5009017845242299296281_0436.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586383, 0.621739478308211, -1.2501828803165005], [-0.3976197158559559, -0.07566485901595878, -2.633910119820667], [1.0879414097563262, -0.08367235852192031, -3.1385123060458495], [1.3412747433710008, -1.2268379978241988, -4.181862535930471], [2.4335785182804317, -0.9652098327338804, -4.878228179559649], [1.479810301087866, -2.396013338547869, -3.5632051489667536], [-0.08698699647460126, -1.4003996926219884, -5.386689524218467], [-1.145899229987119, -2.082013673832742, -4.725979095857028], [-0.23141045828324522, -0.15024159313119034, -6.035999436051609], [0.5934126739618841, -2.390737287525716, -6.404938200328202], [1.350706773188926, 1.0832029916024002, -3.720189049227471], [1.9101736356829002, -0.26833351767859864, -2.113688309236765], [-0.7843720924895164, -1.3445493546988785, -2.4752210448161533], [-1.1581983763610537, 0.4930222734242938, -3.561495300206357], [-0.2560445759534282, 1.8735740976390418, -1.3916701657561936], [-2.0076024771874463, 0.6427130616946733, -1.0543092166280625], [1.5770424436171635, 0.0, 0.0], [2.2927181468939137, 1.391552724358056, 0.0], [2.3410798567223177, 2.0598526928949306, 1.4165023767064728], [1.1453916410070364, 2.033172555891597, 1.9759511228307924], [3.2022179922846443, 1.4422259961598296, 2.1984030628363254], [2.721555559149497, 3.3152059779288, 1.274535884539398], [1.6292323391939711, 2.212255867310578, -0.8090479336198891], [3.54558683009438, 1.2600392214310763, -0.4280914688619843], [1.9974224573334847, -0.6906780683055274, 1.0535722235492997], [1.927718322430896, -0.652993231720624, -1.1102241252095326], [-0.3501493572534746, 0.7132786644586314, 1.0706086973199325], [-0.42668432219275426, -1.2491488329668565, 0.15362238828850308], [0.30845357031669957, -2.19248571212999, -7.307964782882033]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0436', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
