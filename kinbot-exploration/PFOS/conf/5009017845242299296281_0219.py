import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0219'
logfile = 'conf/5009017845242299296281_0219.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863816, 0.7718203945763844, 1.1635336229088489], [-0.346602041513903, 0.2926212946843426, 2.62593772419238], [-0.7363023803695204, -1.1753096862432058, 3.019500458011145], [-0.5090926255515552, -2.175074851556422, 1.8327273118101537], [-0.4446952607212927, -3.4119591218724112, 2.293674406740581], [-1.5014279159495583, -2.0816397791756067, 0.9521362607354191], [1.0867274347625948, -1.814917428111401, 0.913320820343304], [1.4232745804972633, -2.9605502066825755, 0.1408620426916208], [0.9696038499690479, -0.5070465509493643, 0.3831202863898576], [2.052951777282571, -1.7413079281751185, 2.154906415962863], [0.01837074119013288, -1.5602715505038285, 4.04481231065228], [-2.0169579203546344, -1.224069770611329, 3.3629733733434084], [0.9809271994195071, 0.40564959980990223, 2.7245186729056172], [-0.9145558447621683, 1.1282242700023728, 3.4870389756753224], [-2.0119541879597156, 0.6059455318059221, 1.0475948322279935], [-0.4099470658637734, 2.065478369993787, 1.0587059160250747], [1.577042443617166, 0.0, 0.0], [2.2927181468939173, 1.3915527243580537, 0.0], [3.782335574419715, 1.318614735245458, -0.48079377647244637], [4.419592497958774, 0.3492507324843931, 0.14976318527135302], [3.8473329977166375, 1.1104432852891981, -1.7795740852228776], [4.368010411051581, 2.468434085335205, -0.20560554420871863], [2.2938967322202184, 1.8704189044736037, 1.240568989312634], [1.6494649440008835, 2.2352148943143373, -0.8029305726284859], [1.9974224573334824, -0.690678068305529, 1.053572223549297], [1.9277183224308956, -0.6529932317206261, -1.1102241252095328], [-0.3501493572534778, -1.2838136616209423, 0.08241309473864834], [-0.4266843221927563, 0.49153352563555, -1.1586058166012272], [2.7338683532614985, -1.0691945552616584, 2.01145392485856]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0219', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
