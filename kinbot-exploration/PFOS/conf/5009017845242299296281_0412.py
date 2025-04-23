import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0412'
logfile = 'conf/5009017845242299296281_0412.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863837, -1.3935598728846004, 0.08664925740765345], [-0.3976197158559589, -2.2432006455416333, 1.3824827499919214], [-0.749183351766663, -1.5977521800128327, 2.768620777958027], [-0.9314389667275326, -2.683258339050249, 3.885779866230183], [-0.8151262318136152, -2.1246952182318557, 5.077828745612004], [-2.1256098696188777, -3.2593799513859714, 3.7811864484623015], [0.35577392705740013, -4.042446693367131, 3.754047885761596], [0.38701933332052674, -4.746654312881194, 4.98917798458349], [0.1736996060995967, -4.660345229285251, 2.4929344441529544], [1.6380449215560549, -3.1329259859577037, 3.658675933543482], [0.24079307197238695, -0.7855357992989356, 3.1281243020671163], [-1.8795656402123568, -0.909925548699611, 2.670844458158286], [-1.1305489784219331, -3.352053415327129, 1.2479921638695717], [0.8878652548597779, -2.5746029018458136, 1.3820892141482406], [-0.2560445759534282, -2.1420087660532903, -0.9267276815498191], [-2.0076024771874463, -1.2344150958913154, -0.029451230457636982], [1.5770424436171622, 0.0, 0.0], [2.2927181468939146, 1.3915527243580559, 0.0], [2.3410798567223177, 2.0598526928949346, 1.4165023767064728], [1.1453916410070355, 2.033172555891602, 1.975951122830789], [3.2022179922846434, 1.4422259961598285, 2.1984030628363245], [2.721555559149506, 3.3152059779287972, 1.2745358845393995], [1.6292323391939794, 2.212255867310576, -0.8090479336198906], [3.54558683009438, 1.2600392214310714, -0.42809146886198046], [1.9974224573334807, -0.6906780683055236, 1.0535722235493017], [1.927718322430896, -0.6529932317206262, -1.1102241252095302], [-0.35014935725347335, 0.5705349971623056, -1.1530217920585848], [-0.4266843221927601, 0.7576153073313048, 1.0049834283127228], [1.405125201446657, -2.2584903213994383, 3.3165165273762613]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0412', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
