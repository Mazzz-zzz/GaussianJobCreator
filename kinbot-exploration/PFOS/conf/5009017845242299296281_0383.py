import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0383'
logfile = 'conf/5009017845242299296281_0383.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863827, 0.6217394783082115, -1.250182880316503], [-2.2709622836291925, 0.6501421835576484, -1.2334320314121778], [-2.969991788512711, 1.3987418462127592, -0.04493060495238731], [-2.1956407956954154, 2.701984397794216, 0.35697056532747057], [-2.9912731058630286, 3.5033560351048956, 1.0434650438255466], [-1.1300133907845014, 2.399079408117023, 1.0928964051269068], [-1.5840214570309314, 3.6417954560905583, -1.1478019392028769], [-2.622927607728347, 3.651168236684635, -2.1188926501238594], [-0.9345439864168515, 4.80578555502656, -0.6696149084380153], [-0.46936475103302644, 2.6298526165477885, -1.6104624487355241], [-4.201534988011381, 1.7379422988872648, -0.4153843175535555], [-3.0244357979279344, 0.6050007412875054, 1.0168795050389576], [-2.621205471792928, 1.2630031390416208, -2.3676977168360898], [-2.71337798074113, -0.6010536646032272, -1.2662424075749565], [-0.371045161828281, -0.11265730320380662, -2.3156464312139002], [-0.24552532002048857, 1.859884894550717, -1.426765995739983], [1.5770424436171662, 0.0, 0.0], [2.292718146893916, 1.391552724358058, 0.0], [3.782335574419715, 1.3186147352454594, -0.48079377647244775], [4.419592497958777, 0.349250732484398, 0.14976318527135574], [3.847332997716638, 1.1104432852891966, -1.7795740852228787], [4.36801041105158, 2.4684340853352094, -0.20560554420871707], [2.2938967322202157, 1.8704189044736095, 1.2405689893126282], [1.6494649440008797, 2.235214894314337, -0.8029305726284901], [1.9974224573334842, -0.6906780683055245, 1.053572223549302], [1.9277183224308978, -0.6529932317206258, -1.1102241252095282], [-0.35014935725347585, 0.7132786644586373, 1.0706086973199331], [-0.42668432219275615, -1.2491488329668508, 0.1536223882885035], [-0.16762857010986157, 2.096803965174081, -0.8616121305232817]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0383', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
