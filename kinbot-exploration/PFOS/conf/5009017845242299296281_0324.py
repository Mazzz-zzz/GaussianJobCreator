import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0324'
logfile = 'conf/5009017845242299296281_0324.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863819, 0.6217394783082161, -1.2501828803165012], [-0.39761971585595535, -0.0756648590159541, -2.6339101198206674], [-0.7491833517666578, -1.598819837150666, -2.768004365822097], [0.37513150925004024, -2.5102398158058135, -2.164005297273761], [0.8222268376054812, -1.9844600407113069, -1.0369900628863518], [-0.09608992658046997, -3.7311518678641806, -1.926888954299533], [1.8294372022557688, -2.6927750788515765, -3.3359008387581257], [2.942173082714334, -3.1737613002747254, -2.5923269534238123], [1.3352845675176521, -3.3081068248755683, -4.511714083948213], [2.0568917557993647, -1.167436806926696, -3.6555103068320727], [-1.882037352830636, -1.837820181460768, -2.113395933093811], [-0.8963077255402113, -1.9155699118782843, -4.048025307170632], [-1.1305489784219287, 0.5952337900286067, -3.526959494450469], [0.887865254859781, 0.09037708117406636, -2.9207161247297275], [-0.256044575953427, 1.8735740976390451, -1.3916701657561934], [-2.0076024771874446, 0.6427130616946765, -1.054309216628062], [1.577042443617165, 0.0, 0.0], [2.2927181468939164, 1.391552724358056, 0.0], [3.782335574419716, 1.3186147352454591, -0.48079377647244875], [4.419592497958776, 0.3492507324843903, 0.1497631852713494], [3.8473329977166344, 1.1104432852892008, -1.779574085222877], [4.368010411051579, 2.4684340853352054, -0.2056055442087144], [2.293896732220219, 1.8704189044736077, 1.2405689893126342], [1.6494649440008842, 2.2352148943143373, -0.8029305726284879], [1.9974224573334833, -0.6906780683055258, 1.0535722235492986], [1.9277183224308951, -0.652993231720625, -1.110224125209534], [-0.3501493572534751, 0.7132786644586335, 1.0706086973199336], [-0.42668432219275876, -1.2491488329668499, 0.1536223882885011], [1.239325834061783, -0.669322034222557, -3.5161186875600277]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0324', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
