import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0109'
logfile = 'conf/5009017845242299296281_0109.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863853, -1.3935598728845975, 0.08664925740765327], [-2.270962283629195, -1.393254564823212, 0.05367636867327599], [-2.969991788512712, -0.738281968402541, -1.188880669680405], [-3.0567398040787115, 0.8216102792205104, -1.051781739696005], [-3.2568593067482463, 1.3662699044806867, -2.2391317387863072], [-4.051038764929397, 1.1601401368908555, -0.23589746330016825], [-1.4799776082122549, 1.5477397448728933, -0.338886019472651], [-0.3738206562783361, 0.9058213288972065, -0.9608092886381705], [-1.6585634249099084, 2.9517444967194897, -0.2916151561990448], [-1.6211807667533933, 0.9993978241815922, 1.1307824587358108], [-2.265942241823489, -1.0266548422793433, -2.2798191986199297], [-4.2028972120978745, -1.212511974028483, -1.313298995289268], [-2.6212054717929303, -2.681987940783274, 0.09005605494850107], [-2.7133779807411313, -0.7960712600074675, 1.1536489463716018], [-0.3710451618282853, -1.9490799840121085, 1.2553873021032882], [-0.24552532002049143, -2.165558044841978, -0.8973245689258761], [1.5770424436171642, 0.0, 0.0], [2.292718146893919, 1.3915527243580514, 0.0], [3.7823355744197156, 1.3186147352454567, -0.4807937764724427], [4.419592497958776, 0.34925073248439353, 0.14976318527135002], [3.847332997716637, 1.1104432852891968, -1.77957408522288], [4.368010411051582, 2.468434085335203, -0.20560554420871913], [2.293896732220221, 1.8704189044735977, 1.2405689893126346], [1.649464944000881, 2.235214894314338, -0.802930572628487], [1.997422457333482, -0.6906780683055289, 1.0535722235492986], [1.9277183224308945, -0.6529932317206258, -1.110224125209535], [-0.3501493572534748, 0.5705349971623116, -1.1530217920585804], [-0.42668432219275493, 0.7576153073313056, 1.0049834283127257], [-1.2788252973621852, 1.645063722097158, 1.7647151582732417]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0109', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
