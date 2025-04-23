import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0227'
logfile = 'conf/5009017845242299296281_0227.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863803, 0.7718203945763829, 1.163533622908851], [-0.3466020415138996, 0.2926212946843382, 2.6259377241923816], [-0.7363023803695176, -1.1753096862432113, 3.0195004580111435], [-2.2452668735535464, -1.2852675418899087, 3.4324418428267722], [-2.403843994986877, -0.8565477161783767, 4.672423516351344], [-3.003380405095209, -0.5619438307665796, 2.6133162585765817], [-2.864509570633141, -3.054696421585805, 3.3468658361034596], [-1.8669867327177092, -3.902386515343864, 3.9025642260521662], [-4.216284351894492, -3.038084968002563, 3.76851233733842], [-2.855535080987305, -3.223374668817037, 1.780976396438065], [-0.526668395472309, -1.970691804482322, 1.9743151008291824], [0.004125934750529869, -1.5779705333152427, 4.044224788693654], [0.9809271994195121, 0.4056495998098974, 2.7245186729056177], [-0.9145558447621659, 1.1282242700023668, 3.487038975675325], [-2.0119541879597156, 0.6059455318059197, 1.047594832227998], [-0.40994706586377266, 2.0654783699937846, 1.0587059160250827], [1.5770424436171666, 0.0, 0.0], [2.2927181468939177, 1.3915527243580532, 0.0], [2.3410798567223297, 2.0598526928949257, 1.4165023767064695], [1.1453916410070417, 2.033172555891595, 1.9759511228307929], [3.202217992284653, 1.442225996159829, 2.198403062836321], [2.721555559149507, 3.3152059779287977, 1.2745358845394013], [1.6292323391939725, 2.212255867310576, -0.8090479336198861], [3.545586830094381, 1.2600392214310727, -0.4280914688619856], [1.9974224573334816, -0.6906780683055289, 1.0535722235492926], [1.927718322430894, -0.6529932317206233, -1.1102241252095344], [-0.3501493572534775, -1.2838136616209426, 0.08241309473864608], [-0.42668432219275804, 0.4915335256355524, -1.1586058166012265], [-3.5980844797605207, -3.7752803485582525, 1.4981639183772082]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0227', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
