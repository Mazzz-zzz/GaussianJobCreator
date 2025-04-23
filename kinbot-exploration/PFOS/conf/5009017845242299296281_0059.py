import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0059'
logfile = 'conf/5009017845242299296281_0059.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863821, 0.7718203945763772, 1.1635336229088515], [-2.2709622836291916, 0.7431123812655573, 1.1797556627389043], [-2.997023964301902, 1.46568770419547, 2.3681573973696874], [-2.2704853996782592, 2.7960258391239754, 2.7705593355752116], [-1.8789938252023843, 3.4399346516184144, 1.6849177135915088], [-3.0861231913208136, 3.573685520699883, 3.4768740338851982], [-0.7559231095728451, 2.47057410629254, 3.8297314510043576], [0.0545783803479677, 3.6390674808373156, 3.8178328674284883], [-1.2087856252647684, 1.8311883880450632, 5.009332322837855], [-0.07639890773680776, 1.375361648067047, 2.9245799939123613], [-4.239914773111564, 1.7601618341371255, 1.9973488528651564], [-3.0230545829338507, 0.6702887991292378, 3.4298036615283674], [-2.642782557053616, 1.3431159399560006, 0.045483693336772706], [-2.6677893778920323, -0.5232706377056817, 1.1465664052870896], [-0.37104516182828395, 2.061737287215907, 1.0602591291106194], [-0.24552532002048805, 0.305673150291253, 2.3240905646658554], [1.577042443617165, 0.0, 0.0], [2.2927181468939133, 1.3915527243580599, 0.0], [2.3410798567223146, 2.0598526928949354, 1.416502376706473], [1.1453916410070324, 2.0331725558915896, 1.9759511228307978], [3.202217992284636, 1.4422259961598387, 2.1984030628363254], [2.7215555591494893, 3.3152059779288017, 1.2745358845394035], [1.6292323391939627, 2.2122558673105788, -0.8090479336198857], [3.545586830094375, 1.2600392214310863, -0.4280914688619853], [1.997422457333485, -0.6906780683055282, 1.0535722235492953], [1.9277183224308965, -0.6529932317206197, -1.110224125209536], [-0.3501493572534731, -1.2838136616209435, 0.0824130947386427], [-0.4266843221927561, 0.491533525635548, -1.1586058166012245], [0.36090463582346977, 0.709611036508507, 3.47367390272492]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0059', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
