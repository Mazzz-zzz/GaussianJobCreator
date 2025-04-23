import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0056'
logfile = 'conf/5009017845242299296281_0056.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863801, -1.393559872884601, 0.08664925740765349], [-0.3466020415138989, -2.4204394252486665, -1.0595513872112594], [-0.7363023803695169, -2.0273092602547957, -2.527598274606108], [-0.5090926255515532, -0.499650984458944, -2.8000337324856055], [0.6213430942919607, -0.10684176608509977, -2.239393077217541], [-0.45968485817329763, -0.26494129642363523, -4.108227248529034], [-1.8990683064763123, 0.5480559603545136, -2.098648606223011], [-1.4466844030753891, 1.8941772651939206, -2.0229040063416472], [-3.0879999660094115, 0.14642164190718357, -2.754764567660643], [-1.9289234104268815, -0.057258805287837555, -0.6449480397126811], [0.018370741190136732, -2.7227744393129907, -3.373640954864594], [-2.016957920354632, -2.300385488260377, -2.741562204025706], [0.9809271994195139, -2.5623271837262798, -1.0109564779824443], [-0.9145558447621627, -3.583976471922483, -0.7664486088494578], [-2.0119541879597134, -1.2102165034857026, 0.0009668077396031385], [-0.4099470658637678, -1.949605403411485, 1.2594037813693546], [1.5770424436171666, 0.0, 0.0], [2.2927181468939137, 1.3915527243580597, 0.0], [1.6005215470082463, 2.4407219045639, -0.9357086002340298], [1.3760692761371542, 1.914535025410579, -2.1257143081021432], [0.4611709133627442, 2.8519467789832507, -0.4188289776134463], [2.404864907693485, 3.4780800111829153, -1.0689303403306898], [3.5400592233304655, 1.2295174124846087, -0.4315210556927504], [2.308846803952252, 1.8960947387583766, 1.2310220414904676], [1.9974224573334864, -0.6906780683055245, 1.053572223549298], [1.9277183224308962, -0.6529932317206216, -1.1102241252095328], [-0.35014935725347734, 0.5705349971623094, -1.1530217920585795], [-0.42668432219275454, 0.7576153073313002, 1.00498342831273], [-1.5684924715178075, -0.9550596904189905, -0.6472923159006685]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0056', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
