import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0266'
logfile = 'conf/5009017845242299296281_0266.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863816, 0.6217394783082123, -1.2501828803165018], [-0.3976197158559542, -0.07566485901595821, -2.6339101198206674], [1.0879414097563278, -0.08367235852192238, -3.1385123060458486], [1.3412747433710017, -1.2268379978242006, -4.18186253593047], [2.433578518280434, -0.9652098327338812, -4.878228179559648], [1.4798103010878665, -2.396013338547872, -3.5632051489667527], [-0.08698699647460004, -1.4003996926219908, -5.386689524218467], [0.36660723242571114, -2.145547744657815, -6.509826944067679], [-1.2312818756427104, -1.7348476260966852, -4.622431334238047], [-0.21781728104830603, 0.11266548078638995, -5.803894327002465], [1.3507067731889293, 1.0832029916023966, -3.7201890492274705], [1.9101736356829009, -0.2683335176785995, -2.113688309236764], [-0.7843720924895149, -1.3445493546988803, -2.475221044816154], [-1.1581983763610504, 0.4930222734242933, -3.5614953002063583], [-0.2560445759534266, 1.87357409763904, -1.391670165756195], [-2.007602477187446, 0.6427130616946739, -1.0543092166280643], [1.5770424436171635, 0.0, 0.0], [2.2927181468939146, 1.3915527243580548, 0.0], [3.7823355744197142, 1.3186147352454582, -0.4807937764724406], [4.419592497958775, 0.34925073248439453, 0.1497631852713529], [3.847332997716636, 1.110443285289198, -1.7795740852228752], [4.368010411051577, 2.4684340853352063, -0.2056055442087123], [2.293896732220215, 1.870418904473609, 1.2405689893126346], [1.6494649440008788, 2.235214894314338, -0.8029305726284873], [1.9974224573334838, -0.6906780683055282, 1.0535722235493001], [1.927718322430896, -0.6529932317206262, -1.1102241252095317], [-0.3501493572534744, 0.713278664458631, 1.0706086973199318], [-0.42668432219275476, -1.2491488329668525, 0.15362238828850258], [0.13705890072902113, 0.684741437557516, -5.109089043496682]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0266', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
