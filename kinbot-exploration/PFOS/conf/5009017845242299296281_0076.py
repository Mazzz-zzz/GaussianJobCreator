import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0076'
logfile = 'conf/5009017845242299296281_0076.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863848, 0.771820394576383, 1.1635336229088467], [-2.2709622836291934, 0.7431123812655673, 1.1797556627388972], [-3.0203184489305834, 1.3845311280592616, -0.040337232539329373], [-3.1314232385001817, 2.942629924600794, 0.09955209854968268], [-1.9982160731709808, 3.43310750112438, 0.5704822530854788], [-3.396984742281262, 3.4957075334668017, -1.0804827682375686], [-4.510308809776055, 3.4410505216026093, 1.270839774966778], [-4.4625709530279964, 2.583708050158009, 2.404469340345329], [-4.496161508544784, 4.855168652014235, 1.34450481930274], [-5.73434168405041, 3.0226385534204305, 0.37238689507458633], [-2.3454389866601097, 1.108631838906674, -1.1527260342433971], [-4.247573342377321, 0.8876550485557978, -0.12629879149041662], [-2.5961406643712652, -0.55182496350124, 1.2269126406770368], [-2.6897770261666536, 1.3325843308485144, 2.2930463940297656], [-0.37104516182828245, 2.061737287215914, 1.060259129110608], [-0.24552532002049157, 0.30567315029126524, 2.3240905646658527], [1.5770424436171662, 0.0, 0.0], [2.2927181468939164, 1.3915527243580523, 0.0], [3.782335574419716, 1.3186147352454567, -0.4807937764724389], [4.419592497958774, 0.3492507324843902, 0.14976318527135557], [3.847332997716642, 1.1104432852891912, -1.779574085222873], [4.368010411051583, 2.468434085335205, -0.20560554420871696], [2.2938967322202144, 1.8704189044736061, 1.2405689893126353], [1.6494649440008826, 2.2352148943143355, -0.8029305726284865], [1.9974224573334836, -0.6906780683055272, 1.0535722235493024], [1.9277183224308958, -0.6529932317206277, -1.1102241252095264], [-0.3501493572534771, -1.2838136616209452, 0.08241309473865081], [-0.4266843221927556, 0.4915335256355425, -1.158605816601228], [-6.038325866862823, 2.135594173687037, 0.6105218813308377]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0076', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
