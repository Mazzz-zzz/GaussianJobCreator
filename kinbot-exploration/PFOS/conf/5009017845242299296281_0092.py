import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0092'
logfile = 'conf/5009017845242299296281_0092.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863837, 0.7718203945763862, 1.1635336229088493], [-2.2709622836291943, 0.7431123812655629, 1.179755662738897], [-3.0203184489305848, 1.3845311280592572, -0.04033723253932801], [-4.480058160568046, 0.829245790151868, -0.1831774466796777], [-4.45376171327607, -0.35984159707031715, -0.7592691281730998], [-5.056700168542039, 0.7301262577332763, 1.011209154670927], [-5.555811090753358, 1.9566623859997088, -1.2288116005401852], [-4.795842534843593, 2.3780518228087764, -2.354560126803562], [-6.829885671024435, 1.344680560406159, -1.3158114896106898], [-5.672775342070443, 3.1599411412281606, -0.21936510322658176], [-3.086116569243987, 2.7007512384555876, 0.1391809583425341], [-2.3639288743924602, 1.1178871860085935, -1.1621587584935864], [-2.5961406643712652, -0.5518249635012427, 1.226912640677037], [-2.689777026166654, 1.3325843308485115, 2.293046394029765], [-0.37104516182828434, 2.061737287215912, 1.06025912911061], [-0.24552532002049293, 0.3056731502912636, 2.324090564665856], [1.5770424436171668, 0.0, 0.0], [2.2927181468939164, 1.3915527243580534, 0.0], [2.34107985672232, 2.05985269289493, 1.416502376706471], [1.1453916410070368, 2.0331725558915976, 1.9759511228307884], [3.202217992284642, 1.442225996159829, 2.1984030628363294], [2.7215555591495035, 3.3152059779287972, 1.2745358845393995], [1.629232339193976, 2.212255867310573, -0.8090479336198914], [3.54558683009438, 1.2600392214310703, -0.4280914688619817], [1.9974224573334831, -0.6906780683055236, 1.053572223549308], [1.9277183224308985, -0.652993231720628, -1.1102241252095242], [-0.3501493572534741, -1.2838136616209443, 0.08241309473865076], [-0.4266843221927536, 0.4915335256355434, -1.1586058166012274], [-6.549857204053543, 3.563949433201325, -0.27824713862942824]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0092', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
