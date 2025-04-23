import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0402'
logfile = 'conf/5009017845242299296281_0402.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586385, 0.7718203945763881, 1.1635336229088498], [-0.3976197158559556, 2.3188655045575928, 1.2514273698287484], [1.087941409756327, 2.7598675663867454, 1.4967937649483825], [2.1006320219594947, 1.8490216187922408, 0.719360460601157], [2.276257591063596, 0.7136703779757699, 1.3725336095673069], [1.6464483143989896, 1.5907192761912547, -0.5037082619759725], [3.7756073465590045, 2.673466043561305, 0.5285401693191063], [4.111588383345388, 3.289558684377205, 1.765462805809534], [4.610710638282379, 1.7550116729623864, -0.15292944330388322], [3.3619554359701462, 3.794094830612934, -0.49792261534791993], [1.23959077527223, 4.0125250124967184, 1.0762663043113099], [1.3717373297746531, 2.684466857349559, 2.790704275898288], [-0.7843720924895157, 2.8158789821420926, 0.07319662459687586], [-1.158198376361052, 2.837834268725449, 2.2077174635201766], [-0.25604457595342967, 0.26843466841424324, 2.318397847306017], [-2.0076024771874486, 0.5917020341966377, 1.0837604470856956], [1.5770424436171642, 0.0, 0.0], [2.292718146893914, 1.3915527243580537, 0.0], [2.341079856722322, 2.05985269289493, 1.4165023767064735], [1.1453916410070326, 2.033172555891591, 1.9759511228307916], [3.2022179922846403, 1.442225996159828, 2.1984030628363245], [2.7215555591494978, 3.3152059779287955, 1.2745358845394048], [1.6292323391939711, 2.2122558673105757, -0.8090479336198854], [3.545586830094378, 1.2600392214310758, -0.42809146886198135], [1.9974224573334811, -0.6906780683055285, 1.053572223549301], [1.9277183224308918, -0.6529932317206282, -1.1102241252095304], [-0.3501493572534785, -1.2838136616209423, 0.08241309473865063], [-0.42668432219276176, 0.49153352563554564, -1.1586058166012256], [3.161949723043475, 4.618535543792193, -0.03288755006233107]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0402', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
