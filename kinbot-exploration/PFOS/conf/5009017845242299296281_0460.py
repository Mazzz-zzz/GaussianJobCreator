import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0460'
logfile = 'conf/5009017845242299296281_0460.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863816, 0.6217394783082123, -1.2501828803165018], [-0.3976197158559542, -0.07566485901595821, -2.6339101198206674], [1.0879414097563278, -0.08367235852192238, -3.1385123060458486], [1.3412747433710017, -1.2268379978242006, -4.18186253593047], [2.433578518280434, -0.9652098327338812, -4.878228179559648], [1.4798103010878665, -2.396013338547872, -3.5632051489667527], [-0.08698699647460004, -1.4003996926219908, -5.386689524218467], [0.3666072324257096, -2.1455477446578146, -6.509826944067678], [-1.2312818756427102, -1.734847626096685, -4.622431334238046], [-0.21781728104830603, 0.11266548078638995, -5.803894327002465], [1.3507067731889293, 1.0832029916023966, -3.7201890492274705], [1.9101736356829009, -0.2683335176785995, -2.113688309236764], [-0.7843720924895149, -1.3445493546988803, -2.475221044816154], [-1.1581983763610504, 0.4930222734242933, -3.5614953002063583], [-0.2560445759534266, 1.87357409763904, -1.391670165756195], [-2.007602477187446, 0.6427130616946739, -1.0543092166280643], [1.5770424436171635, 0.0, 0.0], [2.2927181468939146, 1.3915527243580548, 0.0], [1.600521547008253, 2.440721904563897, -0.9357086002340322], [1.3760692761371673, 1.9145350254105786, -2.125714308102145], [0.46117091336275595, 2.8519467789832533, -0.41882897761345017], [2.404864907693505, 3.478080011182908, -1.0689303403306858], [3.5400592233304677, 1.2295174124846067, -0.431521055692745], [2.308846803952252, 1.8960947387583729, 1.2310220414904711], [1.9974224573334842, -0.6906780683055285, 1.0535722235493], [1.9277183224308956, -0.652993231720626, -1.1102241252095315], [-0.3501493572534744, 0.713278664458631, 1.0706086973199318], [-0.42668432219275476, -1.2491488329668525, 0.15362238828850258], [-1.1464340732921574, 0.34232821654438084, -5.948422067117977]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0460', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
