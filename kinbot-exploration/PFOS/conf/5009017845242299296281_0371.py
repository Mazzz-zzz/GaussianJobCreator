import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0371'
logfile = 'conf/5009017845242299296281_0371.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586385, 0.6217394783082146, -1.2501828803164983], [-2.270962283629196, 0.6501421835576482, -1.2334320314121707], [-2.9699917885127136, 1.3987418462127532, -0.044930604952376545], [-3.056739804078712, 0.5000645662030728, 1.2374262436633923], [-1.929941607865642, -0.17317330177088175, 1.390979300337153], [-3.275164879501223, 1.2495736306950729, 2.3141930193506854], [-4.462079523721017, -0.7386873768807473, 1.1270757700950764], [-4.465482345408423, -1.2919151885269824, -0.18302285864845463], [-4.424239723399229, -1.5102299556530308, 2.313940113775953], [-5.664988624217497, 0.27019756536760153, 1.2523893729332911], [-2.2659422418234914, 2.4877087631800094, 0.2508004249777487], [-4.2028972120978745, 1.7436062796993261, -0.3934166742568476], [-2.6212054717929334, 1.2630031390416236, -2.3676977168360804], [-2.713377980741131, -0.6010536646032281, -1.2662424075749514], [-0.37104516182828523, -0.11265730320380288, -2.315646431213897], [-0.24552532002049454, 1.8598848945507183, -1.4267659957399748], [1.5770424436171642, 0.0, 0.0], [2.2927181468939133, 1.39155272435806, 0.0], [1.6005215470082366, 2.440721904563899, -0.9357086002340282], [1.3760692761371571, 1.914535025410575, -2.125714308102145], [0.46117091336274485, 2.851946778983253, -0.41882897761344384], [2.4048649076934883, 3.478080011182917, -1.0689303403306851], [3.5400592233304575, 1.2295174124846144, -0.4315210556927483], [2.3088468039522474, 1.8960947387583766, 1.2310220414904682], [1.9974224573334864, -0.6906780683055209, 1.0535722235492986], [1.9277183224308962, -0.652993231720621, -1.1102241252095317], [-0.3501493572534744, 0.7132786644586341, 1.0706086973199365], [-0.42668432219275343, -1.2491488329668503, 0.15362238828850344], [-5.382816026977956, 1.0743907408325106, 1.7102259688916723]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0371', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
