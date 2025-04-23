import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0037'
logfile = 'conf/5009017845242299296281_0037.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863816, 0.7718203945763836, 1.1635336229088498], [-2.270962283629192, 0.7431123812655636, 1.1797556627389005], [-2.9699917885127096, -0.6604598778102136, 1.2338112746327905], [-2.1956407956954127, -1.6601377768739944, 2.1615018464552525], [-1.7814667182728683, -1.0342821878438413, 3.249282411802833], [-2.9823571229144368, -2.6753279933790424, 2.506996108856898], [-0.6934923395346738, -2.3870950746070143, 1.303162499785407], [-1.1396204347114396, -3.3916181622652655, 0.40075788052335587], [0.13508964616168487, -1.29983992498962, 0.9334532839712028], [-0.0343811651830556, -3.0860075756933254, 2.551215545725002], [-4.201534988011381, -0.5092377781085995, 1.7127943399246788], [-3.0244357979279344, -1.1831438545952393, 0.015506258743916439], [-2.6212054717929267, 1.4189848017416529, 2.2776416618875843], [-2.7133779807411282, 1.3971249246106952, 0.11259346120335725], [-0.37104516182828107, 2.061737287215914, 1.0602591291106125], [-0.2455253200204881, 0.30567315029126135, 2.3240905646658554], [1.577042443617165, 0.0, 0.0], [2.2927181468939173, 1.3915527243580539, 0.0], [2.3410798567223265, 2.05985269289493, 1.416502376706473], [1.1453916410070435, 2.0331725558915985, 1.9759511228307907], [3.202217992284646, 1.4422259961598303, 2.1984030628363254], [2.7215555591495106, 3.315205977928799, 1.2745358845393975], [1.629232339193981, 2.212255867310575, -0.809047933619889], [3.5455868300943827, 1.260039221431071, -0.4280914688619861], [1.9974224573334842, -0.6906780683055254, 1.0535722235492981], [1.9277183224308971, -0.6529932317206261, -1.1102241252095308], [-0.3501493572534752, -1.283813661620944, 0.08241309473865073], [-0.4266843221927575, 0.4915335256355464, -1.1586058166012274], [-0.3186719321341855, -2.6526459404501552, 3.3681236431397004]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0037', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
