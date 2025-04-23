import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0012'
logfile = 'conf/5009017845242299296281_0012.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863823, 0.621739478308213, -1.2501828803165018], [-2.270962283629192, 0.6501421835576489, -1.2334320314121774], [-2.997023964301902, 1.3180406141844383, -2.4534014845326206], [-4.4759542789311215, 0.8159291580463162, -2.595649047209702], [-5.168854707529478, 1.6639414173142837, -3.3354234366396462], [-4.503853254513123, -0.38856665038693455, -3.158953170380498], [-5.326480128927349, 0.670886698088345, -0.929146183602041], [-6.72993834376042, 0.6081010601867171, -1.1500398473881641], [-4.598236392479622, -0.2830737496486555, -0.17749646722738566], [-4.967238916287429, 2.0970014225598423, -0.36551050907852134], [-3.0150768858013453, 2.635812688105748, -2.274118452082632], [-2.351095681931021, 1.0276075788966017, -3.5753886090291065], [-2.6427825570536148, -0.6321679360904252, -1.1859143708980975], [-2.6677893778920336, 1.2545909529572579, -0.12011753733594911], [-0.3710451618282827, -0.11265730320380304, -2.3156464312139007], [-0.24552532002048982, 1.8598848945507171, -1.4267659957399765], [1.5770424436171662, 0.0, 0.0], [2.2927181468939155, 1.391552724358056, 0.0], [3.782335574419715, 1.3186147352454607, -0.4807937764724473], [4.419592497958776, 0.34925073248440297, 0.14976318527135085], [3.8473329977166366, 1.1104432852892017, -1.7795740852228783], [4.368010411051577, 2.46843408533521, -0.2056055442087184], [2.2938967322202175, 1.8704189044736066, 1.2405689893126317], [1.6494649440008782, 2.235214894314341, -0.802930572628487], [1.9974224573334838, -0.6906780683055258, 1.0535722235493015], [1.9277183224308947, -0.6529932317206304, -1.1102241252095344], [-0.3501493572534752, 0.713278664458633, 1.0706086973199327], [-0.4266843221927533, -1.2491488329668503, 0.15362238828850006], [-4.169785779052184, 2.4369919177062345, -0.7949721998261113]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0012', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
