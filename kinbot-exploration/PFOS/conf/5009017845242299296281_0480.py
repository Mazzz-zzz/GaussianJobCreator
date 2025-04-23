import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0480'
logfile = 'conf/5009017845242299296281_0480.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863862, 0.6217394783082195, -1.2501828803164967], [-2.2709622836291947, 0.6501421835576613, -1.2334320314121676], [-2.9970239643019054, 1.3180406141844556, -2.4534014845326064], [-3.05161235149174, 2.879165519244264, -2.3137795027818076], [-4.021172178383425, 3.222108939721311, -1.4838123374636436], [-1.8922874678070292, 3.3440017300522937, -1.8568406099918846], [-3.3760113189425525, 3.715776527014534, -3.96196355279244], [-2.166002757627481, 3.710628477991158, -4.7091474571596805], [-4.60814527464748, 3.2108592639035005, -4.443898449079686], [-3.634449254937269, 5.17523602390981, -3.4293000192461793], [-2.332949719823826, 1.0176907712167431, -3.5659602120414684], [-4.241488402883609, 0.8658956272092552, -2.5388589204157697], [-2.6427825570536223, -0.632167936090412, -1.1859143708980893], [-2.6677893778920323, 1.2545909529572672, -0.12011753733593764], [-0.37104516182828834, -0.11265730320379358, -2.315646431213897], [-0.24552532002049088, 1.8598848945507225, -1.42676599573997], [1.5770424436171646, 0.0, 0.0], [2.292718146893918, 1.3915527243580506, 0.0], [1.6005215470082514, 2.4407219045638953, -0.9357086002340266], [1.3760692761371656, 1.9145350254105824, -2.12571430810214], [0.4611709133627677, 2.851946778983256, -0.41882897761344184], [2.4048649076935087, 3.4780800111829078, -1.068930340330682], [3.5400592233304646, 1.2295174124845953, -0.43152105569275073], [2.308846803952262, 1.8960947387583686, 1.231022041490474], [1.9974224573334816, -0.6906780683055321, 1.0535722235492946], [1.9277183224308874, -0.6529932317206271, -1.1102241252095368], [-0.3501493572534701, 0.7132786644586322, 1.0706086973199378], [-0.42668432219276226, -1.2491488329668499, 0.15362238828850225], [-4.3158030237172795, 5.610996959712624, -3.960178381296779]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0480', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
