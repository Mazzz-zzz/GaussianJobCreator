import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0373'
logfile = 'conf/5009017845242299296281_0373.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863831, 0.6217394783082182, -1.250182880316501], [-0.39761971585595723, -0.07566485901595355, -2.6339101198206682], [1.0879414097563247, -0.08367235852191891, -3.138512306045849], [1.3412747433709977, -1.226837997824196, -4.181862535930472], [1.4942413294616155, -2.3804958945597914, -3.555622467014099], [0.3200544491791523, -1.3148826342198838, -5.029363151438025], [2.884934235898715, -0.9110824602207787, -5.201213870219509], [3.9071238216419815, -0.44683395752611815, -4.328268977139902], [3.030909577304207, -2.013109346536125, -6.0784736454216155], [2.3532168904157036, 0.305531680436346, -6.048367532239947], [1.350706773188925, 1.083202991602401, -3.7201890492274705], [1.9101736356828993, -0.2683335176785988, -2.113688309236766], [-0.7843720924895166, -1.3445493546988743, -2.4752210448161533], [-1.1581983763610544, 0.493022273424298, -3.561495300206358], [-0.25604457595342633, 1.8735740976390454, -1.3916701657561923], [-2.007602477187446, 0.6427130616946778, -1.0543092166280625], [1.5770424436171653, 0.0, 0.0], [2.2927181468939146, 1.391552724358057, 0.0], [3.782335574419715, 1.3186147352454585, -0.4807937764724428], [4.419592497958774, 0.34925073248438787, 0.14976318527135013], [3.8473329977166357, 1.1104432852891972, -1.779574085222877], [4.368010411051579, 2.4684340853352027, -0.20560554420871463], [2.2938967322202166, 1.8704189044736084, 1.240568989312635], [1.6494649440008817, 2.2352148943143373, -0.802930572628484], [1.9974224573334867, -0.6906780683055282, 1.0535722235492984], [1.9277183224308942, -0.652993231720624, -1.110224125209533], [-0.3501493572534747, 0.7132786644586312, 1.0706086973199342], [-0.42668432219275654, -1.2491488329668494, 0.1536223882885022], [1.3875907372252267, 0.2809117716503136, -6.102426565360503]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0373', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
