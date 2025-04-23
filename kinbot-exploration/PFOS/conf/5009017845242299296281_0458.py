import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0458'
logfile = 'conf/5009017845242299296281_0458.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586382, 0.6217394783082144, -1.2501828803165003], [-2.270962283629191, 0.6501421835576525, -1.233432031412175], [-2.9970239643019005, 1.3180406141844432, -2.4534014845326193], [-2.2704853996782592, 1.0013618477382615, -3.8067090741066756], [-1.8789938252023834, -0.26078578255258567, -3.8215296524556748], [-3.086123191320812, 1.2242184787530934, -4.833339463005324], [-0.7559231095728453, 2.081357673095724, -4.054445663483564], [-1.1920897005380109, 3.3730607645995607, -4.459035376817032], [0.10363285767357773, 1.8426613177583644, -2.9546527917264584], [-0.14900829291124634, 1.3389267029677108, -5.303842834010374], [-4.239914773111563, 0.8496739297323579, -2.5230192895671464], [-3.02305458293385, 2.6351527013118194, -2.2953889586822815], [-2.6427825570536156, -0.6321679360904207, -1.1859143708980955], [-2.667789377892032, 1.254590952957263, -0.12011753733594736], [-0.3710451618282834, -0.11265730320380345, -2.3156464312138985], [-0.245525320020488, 1.8598848945507176, -1.4267659957399808], [1.5770424436171655, 0.0, 0.0], [2.292718146893918, 1.391552724358054, 0.0], [1.60052154700825, 2.4407219045638993, -0.9357086002340262], [1.376069276137167, 1.914535025410581, -2.1257143081021432], [0.4611709133627575, 2.851946778983253, -0.4188289776134476], [2.4048649076935034, 3.47808001118291, -1.0689303403306827], [3.5400592233304655, 1.2295174124846073, -0.4315210556927446], [2.308846803952251, 1.8960947387583722, 1.2310220414904725], [1.997422457333483, -0.6906780683055276, 1.0535722235493001], [1.9277183224308931, -0.6529932317206264, -1.110224125209531], [-0.35014935725347274, 0.7132786644586345, 1.070608697319933], [-0.42668432219275565, -1.2491488329668485, 0.15362238828850322], [0.8175106075906826, 1.3455180410671217, -5.261888621699798]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0458', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
