import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0246'
logfile = 'conf/5009017845242299296281_0246.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863829, 0.7718203945763857, 1.1635336229088462], [-0.39761971585595646, 2.3188655045575937, 1.251427369828741], [1.087941409756326, 2.7598675663867493, 1.4967937649483751], [1.4730420899455359, 2.685713952893135, 3.015274540831184], [0.9614495666061597, 1.5934621813459509, 3.5554236875418717], [2.795349477668078, 2.6726461892025344, 3.1575153453829055], [0.8248630198382829, 4.1642430621827256, 3.9720811551337936], [-0.5063280479077278, 4.422511558550012, 3.543526919654307], [1.20042992580093, 3.9756363801331513, 5.324387995418162], [1.7738962514906857, 5.263487740628936, 3.362540597060672], [1.8947536492841413, 1.9494330700686122, 0.81764614488947], [1.259810370563856, 4.009816458287682, 1.0862835586981865], [-0.7843720924895171, 2.815878982142088, 0.07319662459686824], [-1.158198376361054, 2.8378342687254494, 2.207717463520167], [-0.25604457595342883, 0.268434668414247, 2.3183978473060134], [-2.0076024771874454, 0.5917020341966354, 1.0837604470856936], [1.5770424436171646, 0.0, 0.0], [2.2927181468939177, 1.3915527243580552, 0.0], [2.3410798567223314, 2.0598526928949292, 1.4165023767064715], [1.14539164100704, 2.0331725558915985, 1.9759511228307893], [3.20221799228465, 1.442225996159828, 2.1984030628363245], [2.7215555591495133, 3.3152059779287937, 1.2745358845394006], [1.629232339193984, 2.212255867310575, -0.8090479336198885], [3.5455868300943854, 1.2600392214310665, -0.4280914688619805], [1.9974224573334853, -0.6906780683055254, 1.0535722235493012], [1.927718322430897, -0.6529932317206293, -1.1102241252095295], [-0.35014935725347696, -1.2838136616209426, 0.08241309473865065], [-0.42668432219275865, 0.4915335256355412, -1.1586058166012287], [2.587602386469016, 4.8586219637296955, 3.030986053598508]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0246', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
