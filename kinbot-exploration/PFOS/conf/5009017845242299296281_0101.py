import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0101'
logfile = 'conf/5009017845242299296281_0101.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863861, 0.6217394783082166, -1.250182880316498], [-2.2709622836291956, 0.6501421835576561, -1.2334320314121698], [-2.9970239643019054, 1.3180406141844492, -2.4534014845326095], [-3.0516123514917433, 2.8791655192442565, -2.3137795027818124], [-4.021172178383429, 3.2221089397213047, -1.4838123374636494], [-1.8922874678070327, 3.3440017300522893, -1.85684060999189], [-3.3760113189425556, 3.7157765270145267, -3.9619635527924464], [-3.8039531646981053, 5.047979006194236, -3.7079294233358615], [-2.2964465318638343, 3.3781628035735354, -4.813969195902862], [-4.621512289586947, 2.8618594756088167, -4.409284054503988], [-2.3329497198238265, 1.0176907712167353, -3.5659602120414715], [-4.241488402883609, 0.8658956272092484, -2.538858920415772], [-2.6427825570536205, -0.6321679360904162, -1.1859143708980908], [-2.667789377892034, 1.2545909529572647, -0.1201175373359402], [-0.37104516182828834, -0.11265730320379719, -2.315646431213897], [-0.2455253200204925, 1.8598848945507214, -1.4267659957399719], [1.577042443617164, 0.0, 0.0], [2.292718146893917, 1.391552724358055, 0.0], [2.341079856722327, 2.059852692894924, 1.4165023767064755], [1.1453916410070453, 2.033172555891593, 1.9759511228307993], [3.2022179922846505, 1.4422259961598263, 2.198403062836322], [2.721555559149511, 3.3152059779287946, 1.2745358845394013], [1.629232339193974, 2.2122558673105743, -0.8090479336198834], [3.5455868300943827, 1.260039221431071, -0.42809146886198635], [1.9974224573334811, -0.6906780683055304, 1.0535722235492957], [1.9277183224308898, -0.6529932317206251, -1.110224125209535], [-0.35014935725347185, 0.7132786644586324, 1.070608697319936], [-0.42668432219276103, -1.2491488329668499, 0.15362238828850225], [-5.442540181202505, 3.3156081122424355, -4.172666720738186]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0101', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
