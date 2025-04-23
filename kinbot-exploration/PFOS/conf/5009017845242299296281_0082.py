import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0082'
logfile = 'conf/5009017845242299296281_0082.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863861, -1.393559872884598, 0.08664925740764336], [-0.39761971585595934, -2.24320064554164, 1.382482749991909], [-1.1233882121466405, -3.6270084497455493, 1.5236046658601465], [-2.5854106473305416, -3.46381850142439, 2.067152822685161], [-3.1611759047412464, -2.4181551748563592, 1.5000685971009284], [-3.2994930518117376, -4.556084848749346, 1.8103638301937368], [-2.6135892170537467, -3.20452056586399, 3.9255483349315394], [-2.3985135550467778, -4.462270636944346, 4.553432262198045], [-1.846409481935056, -2.044084859768806, 4.190465298696755], [-4.132502463634926, -2.817026446654369, 4.078170373796136], [-1.1808982350936317, -4.206246265904724, 0.3277084738360304], [-0.452629532866641, -4.403036426342522, 2.365206219564964], [0.9186702629217368, -2.4703580149136983, 1.3578310147309414], [-0.686339318432844, -1.5018460017288968, 2.4451837914307144], [-0.25604457595343416, -2.142008766053281, -0.9267276815498321], [-2.0076024771874494, -1.2344150958913036, -0.02945123045763835], [1.5770424436171646, 0.0, 0.0], [2.2927181468939173, 1.3915527243580519, 0.0], [1.6005215470082534, 2.4407219045638966, -0.935708600234033], [1.376069276137163, 1.9145350254105784, -2.1257143081021437], [0.4611709133627573, 2.8519467789832573, -0.4188289776134454], [2.4048649076935145, 3.4780800111829056, -1.0689303403306911], [3.5400592233304655, 1.2295174124845953, -0.4315210556927517], [2.3088468039522656, 1.8960947387583704, 1.2310220414904642], [1.9974224573334838, -0.6906780683055355, 1.0535722235492908], [1.9277183224308942, -0.6529932317206283, -1.1102241252095402], [-0.35014935725347546, 0.5705349971623154, -1.1530217920585788], [-0.42668432219275304, 0.7576153073313006, 1.0049834283127266], [-4.655536630702929, -3.5952115987707485, 4.316551404503472]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0082', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
