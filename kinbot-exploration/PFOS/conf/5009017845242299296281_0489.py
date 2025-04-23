import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0489'
logfile = 'conf/5009017845242299296281_0489.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863821, 0.7718203945763845, 1.1635336229088467], [-0.3466020415139028, 0.29262129468434245, 2.6259377241923785], [1.162453549046707, 0.30841610108309503, 3.054984450098301], [1.9243606525365986, -0.9650386511798188, 2.5476426337726172], [1.5424767514304083, -1.2559555411150252, 1.3163693122459377], [3.237845094852695, -0.7569464994859723, 2.5665189880288137], [1.5819802655580883, -2.459696204729642, 3.6295103703023197], [1.9592615238581803, -3.6261934972569145, 2.908834349838176], [2.0501774557344477, -2.143315112978159, 4.927990700357851], [0.009894955985936404, -2.3675313896677634, 3.653924630501017], [1.7462981079693154, 1.3839247474360874, 2.5340254928450237], [1.2535978601361861, 0.3461183998849051, 4.378125052116877], [-1.0117971004258963, 1.1239736000383733, 3.4329051118380076], [-0.8215282967123813, -0.9363359387499617, 2.7884443111439525], [-2.011954187959716, 0.6059455318059169, 1.0475948322279889], [-0.4099470658637771, 2.0654783699937846, 1.0587059160250722], [1.5770424436171684, 0.0, 0.0], [2.2927181468939177, 1.391552724358056, 0.0], [2.3410798567223265, 2.05985269289493, 1.4165023767064748], [1.1453916410070353, 2.03317255589159, 1.975951122830792], [3.202217992284647, 1.4422259961598303, 2.1984030628363262], [2.7215555591495, 3.3152059779287963, 1.2745358845394053], [1.6292323391939774, 2.212255867310576, -0.8090479336198879], [3.545586830094381, 1.2600392214310783, -0.42809146886198013], [1.9974224573334856, -0.6906780683055267, 1.0535722235492986], [1.927718322430896, -0.6529932317206265, -1.1102241252095328], [-0.3501493572534731, -1.2838136616209443, 0.08241309473864848], [-0.42668432219275504, 0.4915335256355466, -1.1586058166012285], [-0.2761471818930415, -1.458985398605793, 3.484553841464536]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0489', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
