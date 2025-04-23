import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0176'
logfile = 'conf/5009017845242299296281_0176.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863806, 0.7718203945763881, 1.1635336229088478], [-0.3466020415139019, 0.29262129468434506, 2.6259377241923803], [-0.9873929842445216, 1.0841115196807807, 3.819465912284933], [-2.4621724959026183, 1.5155496356173241, 3.5056638996836162], [-3.0952492719626608, 0.5339958535318479, 2.887385716691432], [-3.1051333564033308, 1.8185240413419892, 4.629811524283238], [-2.5252703668709384, 3.032185596330106, 2.402261274860259], [-3.8302771896164236, 3.1223350527051075, 1.8443533291263599], [-1.890670372779976, 4.08202950328018, 3.1096984736189675], [-1.5307227375406445, 2.5414473113766984, 1.2839618612489105], [-0.9969507321948188, 0.304121876243912, 4.896738560267488], [-0.27888546472280534, 2.1786256533403408, 4.065404334434963], [-0.7784534903451409, -0.9699349115591062, 2.691387238785716], [0.9731582247379456, 0.30374941440734804, 2.7687950826457883], [-2.011954187959716, 0.6059455318059233, 1.047594832227993], [-0.40994706586377194, 2.065478369993785, 1.0587059160250736], [1.5770424436171664, 0.0, 0.0], [2.292718146893921, 1.3915527243580523, 0.0], [3.782335574419717, 1.318614735245458, -0.48079377647244514], [4.419592497958776, 0.34925073248438765, 0.14976318527135257], [3.84733299771664, 1.1104432852891886, -1.7795740852228783], [4.368010411051585, 2.4684340853352005, -0.20560554420872124], [2.293896732220219, 1.8704189044736026, 1.240568989312631], [1.6494649440008855, 2.2352148943143293, -0.8029305726284899], [1.9974224573334833, -0.6906780683055277, 1.0535722235493015], [1.9277183224308956, -0.6529932317206261, -1.1102241252095315], [-0.3501493572534757, -1.283813661620941, 0.08241309473865055], [-0.4266843221927559, 0.49153352563554537, -1.1586058166012279], [-1.0073517507187166, 3.2836128039730865, 0.9504418791665752]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0176', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
