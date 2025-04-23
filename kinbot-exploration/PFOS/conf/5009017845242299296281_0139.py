import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0139'
logfile = 'conf/5009017845242299296281_0139.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863822, -1.3935598728846004, 0.08664925740765096], [-0.3976197158559557, -2.243200645541639, 1.3824827499919166], [-0.7491833517666574, -1.5977521800128431, 2.7686207779580227], [0.3751315092500395, -0.6189636534602693, 3.255934098715885], [-0.10815375808507674, 0.19767795751880263, 4.175735040589357], [1.3965384052652348, -1.3030408487940413, 3.7635425936703064], [1.0485429431935351, 0.42508470237727286, 1.8494692387133023], [-0.048222455147847995, 0.869431883668178, 1.0607409893841113], [2.0248003447344605, 1.285943845961471, 2.407322079964507], [1.806813834672172, -0.7111493844981632, 1.0655432315113846], [-1.882037352830636, -0.91134447558358, 2.6482969312796554], [-0.8963077255402094, -2.5479077952329394, 3.682944859997024], [-1.1305489784219278, -3.3520534153271333, 1.2479921638695664], [0.8878652548597809, -2.5746029018458185, 1.3820892141482315], [-0.2560445759534267, -2.142008766053288, -0.926727681549826], [-2.0076024771874446, -1.2344150958913145, -0.02945123045763531], [1.5770424436171646, 0.0, 0.0], [2.292718146893914, 1.391552724358056, 0.0], [3.7823355744197142, 1.3186147352454618, -0.48079377647244187], [4.419592497958776, 0.3492507324843962, 0.14976318527135418], [3.8473329977166366, 1.1104432852892014, -1.7795740852228767], [4.368010411051579, 2.4684340853352076, -0.20560554420871197], [2.293896732220217, 1.8704189044736044, 1.2405689893126328], [1.6494649440008784, 2.2352148943143386, -0.802930572628483], [1.9974224573334842, -0.6906780683055277, 1.0535722235492966], [1.9277183224308962, -0.6529932317206237, -1.110224125209532], [-0.350149357253478, 0.5705349971623094, -1.153021792058581], [-0.4266843221927606, 0.7576153073313022, 1.0049834283127268], [1.2485653089482147, -1.0605853313087523, 0.3568715943622641]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0139', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
