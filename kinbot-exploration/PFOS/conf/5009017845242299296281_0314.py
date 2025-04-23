import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0314'
logfile = 'conf/5009017845242299296281_0314.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863883, 0.7718203945763867, 1.1635336229088433], [-0.3976197158559593, 2.318865504557597, 1.2514273698287322], [-1.1233882121466403, 3.1329845708321886, 2.379279124290352], [-2.585410647330543, 3.522116108662275, 1.9661784049894433], [-3.1611759047412464, 2.5081750999368824, 1.344149513167917], [-3.2994930518117367, 3.845863491414979, 3.0405033057174045], [-2.6135892170537485, 5.001884864766419, 0.8124220495220773], [-1.5773396053767719, 4.837209363086736, -0.14752572471142328], [-3.9699523510745336, 5.240341942942602, 0.48261294991293996], [-2.173316102013871, 6.105120084498622, 1.8466242850275538], [-1.1808982350936326, 2.3869269963298305, 3.4788618839288854], [-0.452629532866642, 4.249846884503502, 2.6305382892183458], [0.9186702629217364, 2.411095160260261, 1.4604772899922631], [-0.6863393184328458, 2.8685142811653974, 0.07804489435392167], [-0.25604457595343455, 0.2684346684142552, 2.318397847306012], [-2.00760247718745, 0.5917020341966388, 1.083760447085686], [1.5770424436171646, 0.0, 0.0], [2.2927181468939146, 1.3915527243580526, 0.0], [1.6005215470082519, 2.4407219045638984, -0.9357086002340286], [1.376069276137168, 1.9145350254105749, -2.125714308102143], [0.4611709133627475, 2.851946778983247, -0.41882897761345517], [2.404864907693498, 3.47808001118291, -1.0689303403306807], [3.5400592233304664, 1.2295174124846093, -0.4315210556927379], [2.3088468039522425, 1.8960947387583746, 1.2310220414904713], [1.9974224573334802, -0.690678068305524, 1.0535722235493075], [1.927718322430898, -0.6529932317206335, -1.110224125209525], [-0.35014935725347524, -1.2838136616209448, 0.08241309473865537], [-0.42668432219275354, 0.4915335256355338, -1.158605816601232], [-2.3921091252147035, 5.827869345286887, 2.7473040434884575]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0314', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
