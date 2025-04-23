import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0487'
logfile = 'conf/5009017845242299296281_0487.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863862, 0.7718203945763923, 1.1635336229088384], [-0.39761971585595773, 2.3188655045576003, 1.25142736982873], [-1.123388212146639, 3.132984570832193, 2.3792791242903473], [-2.5854106473305394, 3.5221161086622828, 1.9661784049894373], [-3.290357209347752, 3.8253563239326023, 3.0420955298740266], [-2.5711488993453666, 4.562809768561843, 1.1381365284262568], [-3.456801175356165, 2.1167294688457625, 1.0789184253851205], [-2.9711671724751763, 2.064936044637663, -0.25671178701072894], [-3.4619816816099953, 1.0117689482496615, 1.9645776047125563], [-4.910270578146565, 2.7231666842170696, 1.064662592409996], [-1.1808982350936326, 2.3869269963298367, 3.478861883928881], [-0.4526295328666377, 4.249846884503507, 2.630538289218343], [0.9186702629217396, 2.411095160260262, 1.4604772899922622], [-0.686339318432839, 2.8685142811654, 0.0780448943539192], [-0.2560445759534377, 0.2684346684142572, 2.3183978473060116], [-2.00760247718745, 0.5917020341966445, 1.0837604470856828], [1.5770424436171655, 0.0, 0.0], [2.2927181468939195, 1.3915527243580543, 0.0], [3.782335574419719, 1.3186147352454551, -0.4807937764724365], [4.419592497958775, 0.34925073248438787, 0.14976318527136345], [3.847332997716647, 1.11044328528918, -1.7795740852228685], [4.368010411051582, 2.468434085335201, -0.20560554420871846], [2.2938967322202144, 1.8704189044736128, 1.2405689893126284], [1.6494649440008882, 2.23521489431433, -0.8029305726284962], [1.9974224573334785, -0.6906780683055238, 1.05357222354931], [1.9277183224308985, -0.6529932317206366, -1.110224125209522], [-0.3501493572534777, -1.2838136616209437, 0.08241309473865416], [-0.4266843221927506, 0.4915335256355372, -1.1586058166012345], [-5.015437645029352, 3.3608086173331007, 1.784603589903069]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0487', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
