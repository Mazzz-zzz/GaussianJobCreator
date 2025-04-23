import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0044'
logfile = 'conf/5009017845242299296281_0044.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863808, -1.393559872884599, 0.08664925740765088], [-2.2709622836291916, -1.3932545648232166, 0.05367636867327865], [-2.997023964301899, -2.783728318379926, 0.08524408716293909], [-4.475954278931121, -2.6558625932156357, 0.591209145048289], [-4.4929097750280675, -2.5598963502439855, 1.909163975245065], [-5.055401694144313, -1.5823291229684167, 0.06142935889679178], [-5.510548804511903, -4.144001785101999, 0.10467112564480108], [-4.736207008264467, -5.315917935989797, 0.3270184343368993], [-6.8059689138757, -3.9361359892358667, 0.6376207694319428], [-5.583541308689301, -3.87341972470764, -1.4451675991055069], [-3.0150768858013444, -3.287350694771377, -1.1456215214756182], [-2.351095681931018, -3.6101811532690204, 0.8977600360686682], [-2.6427825570536134, -0.7109480038655857, 1.1404306775613344], [-2.6677893778920354, -0.7313203152515859, -1.0264488679511359], [-0.371045161828279, -1.9490799840121134, 1.2553873021032862], [-0.24552532002049007, -2.1655580448419793, -0.8973245689258789], [1.5770424436171666, 0.0, 0.0], [2.2927181468939155, 1.391552724358058, 0.0], [1.6005215470082472, 2.440721904563901, -0.935708600234026], [1.3760692761371565, 1.9145350254105802, -2.1257143081021406], [0.4611709133627544, 2.851946778983253, -0.41882897761343996], [2.404864907693495, 3.478080011182911, -1.0689303403306816], [3.54005922333046, 1.2295174124846104, -0.4315210556927496], [2.308846803952251, 1.8960947387583755, 1.231022041490469], [1.9974224573334867, -0.6906780683055254, 1.0535722235492966], [1.9277183224308936, -0.652993231720624, -1.1102241252095388], [-0.35014935725347746, 0.5705349971623116, -1.1530217920585761], [-0.4266843221927545, 0.7576153073313003, 1.004983428312728], [-5.456557680368095, -2.931983479011913, -1.6282941676987204]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0044', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
