import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0254'
logfile = 'conf/5009017845242299296281_0254.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863837, 0.77182039457638, 1.1635336229088535], [-0.3976197158559589, 2.318865504557587, 1.2514273698287544], [-0.749183351766663, 3.1965720171635037, -0.0006164121359210983], [-0.9314389667275326, 4.706813247194559, 0.38087995341887876], [-0.8151262318136152, 5.459876298882796, -0.6988743385178752], [-2.1256098696188777, 4.904293496506794, 0.9321126142548074], [0.3557739270574001, 5.272324182776364, 1.62383758701956], [0.38701933332052674, 6.694082035091941, 1.6161342256463604], [0.1736996060995967, 4.48911717324831, 2.789510136890183], [1.6380449215560549, 4.734969295642251, 0.8838555252440579], [0.24079307197238695, 3.101803011435061, -0.8837681932585538], [-1.8795656402123568, 2.767981924671768, -0.5474035883527751], [-1.130548978421933, 2.7568196252985184, 2.2789673305809055], [0.8878652548597777, 2.4842258206717482, 1.5386269105814974], [-0.2560445759534282, 0.26843466841423713, 2.3183978473060183], [-2.0076024771874463, 0.5917020341966307, 1.0837604470857034], [1.5770424436171622, 0.0, 0.0], [2.2927181468939115, 1.391552724358057, 0.0], [1.6005215470082497, 2.4407219045639055, -0.9357086002340256], [1.3760692761371576, 1.914535025410581, -2.1257143081021397], [0.4611709133627502, 2.8519467789832547, -0.41882897761344384], [2.4048649076934936, 3.4780800111829144, -1.0689303403306825], [3.540059223330461, 1.229517412484608, -0.4315210556927502], [2.3088468039522505, 1.8960947387583742, 1.2310220414904707], [1.997422457333484, -0.6906780683055261, 1.0535722235492955], [1.927718322430891, -0.6529932317206223, -1.110224125209535], [-0.3501493572534754, -1.2838136616209457, 0.08241309473865084], [-0.426684322192761, 0.49153352563554603, -1.158605816601223], [1.405125201446657, 4.001432725478509, 0.29765172884508156]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0254', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
