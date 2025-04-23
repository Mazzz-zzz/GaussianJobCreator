import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0345'
logfile = 'conf/5009017845242299296281_0345.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863801, 0.7718203945763835, 1.163533622908852], [-0.3466020415138989, 0.29262129468433806, 2.6259377241923807], [-0.7363023803695169, -1.1753096862432113, 3.0195004580111435], [-0.5090926255515532, -2.1750748515564253, 1.83272731181015], [0.6213430942919607, -1.885950410886849, 1.21222422222366], [-0.45968485817329763, -3.4253585135337774, 2.283559517478961], [-1.8990683064763123, -2.0915109867831907, 0.5746939187490164], [-1.4466844030753894, -2.698974891506144, -0.6289536277580484], [-3.087999966009413, -2.4589069179929646, 1.2505774222748698], [-1.9289234104268815, -0.5299119838682383, 0.37206159982595416], [0.018370741190136732, -1.5602715505038391, 4.044812310652276], [-2.016957920354632, -1.2240697706113324, 3.3629733733434075], [0.9809271994195139, 0.40564959980989734, 2.7245186729056172], [-0.9145558447621627, 1.128224270002364, 3.4870389756753264], [-2.0119541879597134, 0.6059455318059214, 1.0475948322279964], [-0.4099470658637678, 2.065478369993785, 1.0587059160250805], [1.5770424436171666, 0.0, 0.0], [2.292718146893921, 1.3915527243580514, 0.0], [3.782335574419717, 1.3186147352454576, -0.48079377647244514], [4.419592497958776, 0.3492507324843892, 0.14976318527135046], [3.8473329977166406, 1.1104432852892003, -1.7795740852228756], [4.368010411051582, 2.4684340853352023, -0.20560554420871946], [2.293896732220222, 1.8704189044736013, 1.240568989312636], [1.6494649440008862, 2.2352148943143377, -0.8029305726284814], [1.9974224573334824, -0.6906780683055287, 1.0535722235492955], [1.9277183224308945, -0.6529932317206233, -1.1102241252095357], [-0.35014935725347757, -1.2838136616209403, 0.08241309473864708], [-0.42668432219275454, 0.49153352563555275, -1.1586058166012252], [-1.4002647520976845, -0.28091617252734785, -0.3989659584391631]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0345', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
