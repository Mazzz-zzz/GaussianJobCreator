import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0106'
logfile = 'conf/5009017845242299296281_0106.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863837, -1.393559872884597, 0.08664925740765075], [-2.270962283629194, -1.3932545648232142, 0.05367636867328103], [-2.9699917885127154, -0.7382819684025402, -1.188880669680398], [-3.0567398040787155, 0.8216102792205107, -1.0517817396959925], [-4.057927071959131, 1.147154978605423, -0.2529409255073734], [-1.9208515428855029, 1.3069791188858964, -0.5587553832825138], [-3.3438782542772696, 1.6524617020879784, -2.709762666803322], [-4.321509485593893, 0.9024001810351716, -3.419745693361205], [-3.437187055756235, 3.043576208758483, -2.4618702132161774], [-1.9318779018872667, 1.3759659063727214, -3.3503465485949975], [-2.265942241823496, -1.0266548422793405, -2.2798191986199243], [-4.202897212097877, -1.2125119740284827, -1.3132989952892582], [-2.6212054717929294, -2.681987940783275, 0.09005605494850588], [-2.713377980741128, -0.7960712600074716, 1.153648946371608], [-0.37104516182828085, -1.9490799840121122, 1.2553873021032855], [-0.24552532002049343, -2.1655580448419767, -0.8973245689258768], [1.5770424436171644, 0.0, 0.0], [2.2927181468939164, 1.391552724358053, 0.0], [1.6005215470082557, 2.440721904563902, -0.9357086002340227], [1.3760692761371656, 1.9145350254105828, -2.125714308102139], [0.46117091336275484, 2.8519467789832573, -0.4188289776134384], [2.404864907693502, 3.4780800111829118, -1.0689303403306771], [3.540059223330463, 1.229517412484604, -0.4315210556927499], [2.308846803952259, 1.896094738758369, 1.2310220414904733], [1.9974224573334876, -0.6906780683055289, 1.0535722235492897], [1.927718322430891, -0.6529932317206212, -1.1102241252095384], [-0.35014935725347857, 0.5705349971623135, -1.1530217920585748], [-0.4266843221927536, 0.7576153073313019, 1.0049834283127301], [-1.9617652723712191, 0.5794803910424058, -3.8986894665309877]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0106', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
