import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0335'
logfile = 'conf/5009017845242299296281_0335.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863822, -1.3935598728846004, 0.08664925740765096], [-0.3976197158559557, -2.243200645541639, 1.3824827499919166], [-0.7491833517666574, -1.5977521800128431, 2.7686207779580227], [0.3751315092500395, -0.6189636534602693, 3.255934098715885], [-0.10815375808507674, 0.19767795751880263, 4.175735040589357], [1.3965384052652348, -1.3030408487940413, 3.7635425936703064], [1.0485429431935351, 0.42508470237727286, 1.8494692387133023], [-0.0482224551478483, 0.869431883668178, 1.0607409893841113], [2.0248003447344605, 1.2859438459614707, 2.4073220799645068], [1.806813834672172, -0.7111493844981632, 1.0655432315113846], [-1.882037352830636, -0.91134447558358, 2.6482969312796554], [-0.8963077255402094, -2.5479077952329394, 3.682944859997024], [-1.1305489784219278, -3.3520534153271333, 1.2479921638695664], [0.8878652548597809, -2.5746029018458185, 1.3820892141482315], [-0.2560445759534267, -2.142008766053288, -0.926727681549826], [-2.0076024771874446, -1.2344150958913145, -0.02945123045763531], [1.5770424436171646, 0.0, 0.0], [2.292718146893914, 1.391552724358056, 0.0], [2.341079856722317, 2.059852692894931, 1.4165023767064735], [1.1453916410070317, 2.0331725558915967, 1.9759511228307918], [3.202217992284643, 1.4422259961598307, 2.1984030628363267], [2.7215555591494924, 3.315205977928801, 1.2745358845394026], [1.6292323391939707, 2.212255867310577, -0.8090479336198834], [3.5455868300943765, 1.2600392214310787, -0.4280914688619857], [1.9974224573334842, -0.6906780683055278, 1.0535722235492964], [1.927718322430896, -0.6529932317206238, -1.1102241252095317], [-0.350149357253478, 0.5705349971623094, -1.153021792058581], [-0.4266843221927606, 0.7576153073313022, 1.0049834283127268], [2.0364953443706857, -1.4362105044386533, 1.6634522238749834]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0335', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
