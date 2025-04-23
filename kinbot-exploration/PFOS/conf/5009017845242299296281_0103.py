import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0103'
logfile = 'conf/5009017845242299296281_0103.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863831, 0.7718203945763837, 1.1635336229088455], [-0.3976197158559577, 2.3188655045575928, 1.2514273698287433], [-0.74918335176666, 3.1965720171635033, -0.0006164121359367134], [-0.9314389667275313, 4.70681324719456, 0.3808799534188578], [-2.1269824923960123, 4.894694865930261, 0.9118019049808777], [0.003036253824339577, 5.08081870789914, 1.250233746814102], [-0.7791576068935264, 5.827600918475933, -1.116554090025386], [0.6003697055768612, 5.9568696807394, -1.4369054425449233], [-1.7728499369113317, 5.41872935264048, -2.038918642446523], [-1.2703030198852923, 7.166951975057827, -0.44912801895899623], [0.24079307197239408, 3.101803011435057, -0.8837681932585663], [-1.8795656402123533, 2.7679819246717674, -0.5474035883527918], [-1.1305489784219354, 2.7568196252985278, 2.2789673305808926], [0.8878652548597774, 2.4842258206717545, 1.5386269105814903], [-0.25604457595342883, 0.268434668414247, 2.3183978473060134], [-2.007602477187447, 0.5917020341966331, 1.0837604470856919], [1.577042443617165, 0.0, 0.0], [2.292718146893915, 1.3915527243580572, 0.0], [3.782335574419717, 1.3186147352454696, -0.4807937764724428], [4.419592497958775, 0.34925073248440075, 0.14976318527135724], [3.8473329977166406, 1.1104432852891994, -1.7795740852228756], [4.368010411051578, 2.468434085335212, -0.20560554420871946], [2.2938967322202126, 1.8704189044736081, 1.2405689893126326], [1.6494649440008802, 2.235214894314339, -0.8029305726284898], [1.9974224573334818, -0.6906780683055255, 1.0535722235493017], [1.9277183224308994, -0.6529932317206267, -1.1102241252095295], [-0.3501493572534727, -1.283813661620946, 0.08241309473865085], [-0.4266843221927559, 0.49153352563554176, -1.1586058166012287], [-1.8331128165883035, 6.975413197425342, 0.31410223207548976]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0103', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
