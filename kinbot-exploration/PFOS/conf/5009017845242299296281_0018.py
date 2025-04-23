import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0018'
logfile = 'conf/5009017845242299296281_0018.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863842, 0.7718203945763854, 1.1635336229088438], [-0.346602041513907, 0.2926212946843472, 2.6259377241923785], [-0.9873929842445301, 1.0841115196807833, 3.8194659122849295], [-0.18477695124075494, 2.390826822806657, 4.1479286638168755], [0.18832441293676877, 2.9824134714161876, 3.026624712021103], [-0.9332368412513804, 3.2267902686431946, 4.861982548914765], [1.3577185656456103, 2.027649090058602, 5.153151448047535], [1.9850148400569665, 0.8752696424410403, 4.604560525202218], [2.014089937717091, 3.265176117994472, 5.360661083248477], [0.658322754877211, 1.6261622674354417, 6.505998986001677], [-2.2305068985543404, 1.4246952716906192, 3.491346945411261], [-1.0036314587762358, 0.3235218506680432, 4.906512125113517], [-0.7784534903451404, -0.969934911559108, 2.6913872387857136], [0.9731582247379402, 0.3037494144073553, 2.7687950826457866], [-2.011954187959719, 0.6059455318059167, 1.0475948322279884], [-0.40994706586378044, 2.065478369993784, 1.0587059160250685], [1.5770424436171637, 0.0, 0.0], [2.2927181468939146, 1.3915527243580605, 0.0], [2.341079856722319, 2.059852692894936, 1.416502376706473], [1.1453916410070297, 2.0331725558916003, 1.9759511228307909], [3.202217992284643, 1.4422259961598338, 2.1984030628363294], [2.721555559149499, 3.3152059779288017, 1.274535884539403], [1.6292323391939756, 2.2122558673105788, -0.8090479336198879], [3.545586830094379, 1.2600392214310723, -0.42809146886197724], [1.997422457333483, -0.6906780683055198, 1.0535722235493], [1.9277183224308985, -0.6529932317206236, -1.1102241252095286], [-0.3501493572534714, -1.2838136616209426, 0.08241309473864836], [-0.4266843221927569, 0.4915335256355429, -1.1586058166012314], [1.185364229486392, 1.9276642208050832, 7.259183055703147]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0018', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
