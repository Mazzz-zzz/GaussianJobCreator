import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0300'
logfile = 'conf/5009017845242299296281_0300.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863803, -1.3935598728846015, 0.08664925740765352], [-2.2709622836291934, -1.3932545648232195, 0.05367636867327256], [-3.0203184489305865, -0.6573324959322197, 1.2192077454993024], [-2.2823447039971096, 0.6594581441991478, 1.6446984279510868], [-3.113358927133044, 1.4464692125745295, 2.305403265437101], [-1.2352156515326078, 0.37551323269565634, 2.4139812597295442], [-1.6404353738143176, 1.6101510970989743, 0.15950529160622698], [-0.4619365315164053, 0.9647019943859115, -0.3062918126348624], [-2.760091515584602, 1.9027399695106941, -0.6566556285467378], [-1.2298799844268766, 2.9383442180890627, 0.8997011679852719], [-4.2454188440806915, -0.34008637241846795, 0.8099700744271588], [-3.0940580938002635, -1.452014544042186, 2.279147700720935], [-2.596140664371262, -0.786625033299948, -1.0913507571730212], [-2.6897770261666483, -2.6521285947103426, 0.007528686185008044], [-0.3710451618282825, -1.9490799840121102, 1.2553873021032895], [-0.2455253200204851, -2.165558044841984, -0.8973245689258742], [1.5770424436171646, 0.0, 0.0], [2.292718146893914, 1.3915527243580574, 0.0], [3.782335574419715, 1.3186147352454647, -0.4807937764724365], [4.419592497958775, 0.34925073248439586, 0.1497631852713559], [3.8473329977166406, 1.1104432852892017, -1.7795740852228725], [4.368010411051575, 2.4684340853352134, -0.20560554420871074], [2.2938967322202144, 1.8704189044736086, 1.2405689893126335], [1.649464944000881, 2.235214894314338, -0.8029305726284879], [1.9974224573334816, -0.690678068305525, 1.0535722235493024], [1.9277183224308974, -0.652993231720624, -1.1102241252095286], [-0.35014935725347524, 0.5705349971623062, -1.1530217920585835], [-0.4266843221927603, 0.7576153073313034, 1.004983428312721], [-1.7392080181026446, 3.039361384995512, 1.7159994423912457]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0300', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
