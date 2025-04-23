import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0146'
logfile = 'conf/5009017845242299296281_0146.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863813, 0.7718203945763852, 1.16353362290885], [-0.34660204151389995, 0.2926212946843426, 2.62593772419238], [-0.7363023803695175, -1.175309686243206, 3.0195004580111453], [-0.5090926255515525, -2.175074851556422, 1.8327273118101537], [0.6213430942919618, -1.8859504108868477, 1.2122242222236614], [-0.4596848581732949, -3.4253585135337716, 2.2835595174789702], [-1.8990683064763114, -2.0915109867831903, 0.5746939187490203], [-3.018501756452613, -2.7927897717250154, 1.1014840331918252], [-1.94865260897711, -0.7565540325474923, 0.1048129775189399], [-1.2380964843927418, -2.96976167913623, -0.5532841845854125], [0.018370741190135775, -1.5602715505038358, 4.044812310652278], [-2.016957920354632, -1.2240697706113293, 3.362973373343409], [0.9809271994195116, 0.40564959980990234, 2.724518672905618], [-0.9145558447621639, 1.1282242700023726, 3.487038975675322], [-2.0119541879597134, 0.605945531805922, 1.0475948322279953], [-0.4099470658637693, 2.0654783699937864, 1.0587059160250767], [1.5770424436171675, 0.0, 0.0], [2.292718146893919, 1.3915527243580543, 0.0], [3.782335574419717, 1.3186147352454591, -0.48079377647244637], [4.419592497958778, 0.3492507324843882, 0.14976318527134924], [3.847332997716642, 1.1104432852891981, -1.7795740852228783], [4.368010411051579, 2.468434085335205, -0.20560554420872068], [2.2938967322202206, 1.870418904473605, 1.2405689893126333], [1.649464944000885, 2.2352148943143377, -0.8029305726284846], [1.997422457333484, -0.6906780683055277, 1.0535722235492961], [1.9277183224308945, -0.6529932317206232, -1.1102241252095344], [-0.35014935725347607, -1.2838136616209415, 0.08241309473864829], [-0.4266843221927541, 0.49153352563555025, -1.158605816601227], [-1.4776687565480118, -2.630500509686696, -1.4270515289134655]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0146', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
