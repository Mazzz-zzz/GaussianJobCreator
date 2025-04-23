import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0309'
logfile = 'conf/5009017845242299296281_0309.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863809, -1.3935598728845988, 0.08664925740765708], [-2.270962283629191, -1.3932545648232177, 0.05367636867327993], [-2.9970239643018997, -2.783728318379926, 0.08524408716294404], [-4.4759542789311215, -2.655862593215633, 0.5912091450482932], [-5.036363599701665, -1.5773343835482123, 0.07234477238229249], [-5.178548597123876, -3.730553107169243, 0.24436398677965085], [-4.5636267180482815, -2.5129142345826603, 2.4603112399854252], [-3.5410856494537764, -1.6212383872619194, 2.886638525239893], [-5.9315139537235275, -2.3710403132154703, 2.798117504818753], [-4.133463285573643, -3.984341420456059, 2.821426301763957], [-3.015076885801341, -3.287350694771383, -1.145621521475609], [-2.3510956819310196, -3.6101811532690156, 0.8977600360686806], [-2.6427825570536125, -0.7109480038655842, 1.1404306775613366], [-2.667789377892033, -0.7313203152515897, -1.0264488679511392], [-0.3710451618282821, -1.949079984012109, 1.2553873021032933], [-0.24552532002048846, -2.165558044841982, -0.8973245689258721], [1.5770424436171666, 0.0, 0.0], [2.2927181468939164, 1.3915527243580572, 0.0], [2.3410798567223225, 2.0598526928949275, 1.4165023767064746], [1.1453916410070366, 2.0331725558915936, 1.9759511228307944], [3.202217992284645, 1.442225996159825, 2.1984030628363276], [2.721555559149503, 3.3152059779287946, 1.2745358845394037], [1.6292323391939747, 2.212255867310576, -0.8090479336198813], [3.54558683009438, 1.2600392214310723, -0.4280914688619733], [1.9974224573334851, -0.6906780683055234, 1.0535722235493015], [1.9277183224308942, -0.6529932317206295, -1.1102241252095326], [-0.3501493572534734, 0.5705349971623054, -1.1530217920585804], [-0.426684322192755, 0.7576153073313038, 1.004983428312727], [-4.6120449169254965, -4.291917948416297, 3.603934249504121]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0309', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
