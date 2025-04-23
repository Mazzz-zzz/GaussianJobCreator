import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0164'
logfile = 'conf/5009017845242299296281_0164.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863852, 0.7718203945763857, 1.1635336229088462], [-0.3976197158559584, 2.318865504557598, 1.2514273698287353], [-1.123388212146639, 3.132984570832187, 2.3792791242903557], [-1.2235971086413808, 2.3151670007983443, 3.7137144948532206], [-2.214848248727895, 1.444821069257295, 3.6335828543850393], [-0.08476421785157753, 1.6691858167858273, 3.9475263793033895], [-1.546098718177794, 3.4271873548415677, 5.190558012344154], [-2.0148297746218793, 2.618454745935847, 6.262320006364917], [-0.4447116237159659, 4.3120314719147625, 5.287274754854862], [-2.760362629999281, 4.239824636692158, 4.602572952667068], [-0.43435958530525626, 4.245343424397677, 2.6177930333052006], [-2.35331956740016, 3.4429499950288975, 1.9898948685020796], [0.9186702629217383, 2.4110951602602606, 1.4604772899922658], [-0.6863393184328432, 2.8685142811654005, 0.07804489435392688], [-0.25604457595343316, 0.26843466841425334, 2.318397847306014], [-2.00760247718745, 0.59170203419664, 1.0837604470856883], [1.5770424436171635, 0.0, 0.0], [2.292718146893915, 1.3915527243580539, 0.0], [2.3410798567223225, 2.0598526928949292, 1.4165023767064746], [1.1453916410070366, 2.0331725558915985, 1.9759511228307898], [3.2022179922846394, 1.4422259961598234, 2.1984030628363307], [2.721555559149509, 3.3152059779287937, 1.2745358845394037], [1.6292323391939838, 2.212255867310576, -0.8090479336198899], [3.545586830094384, 1.2600392214310683, -0.4280914688619771], [1.9974224573334831, -0.690678068305522, 1.0535722235493097], [1.9277183224308971, -0.6529932317206324, -1.1102241252095253], [-0.3501493572534775, -1.2838136616209472, 0.08241309473865668], [-0.42668432219275687, 0.49153352563553726, -1.1586058166012314], [-2.7315131263772168, 4.2345686847259785, 3.635565988987027]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0164', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
