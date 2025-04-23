import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0148'
logfile = 'conf/5009017845242299296281_0148.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863847, 0.7718203945763836, 1.1635336229088478], [-2.270962283629194, 0.7431123812655617, 1.1797556627388976], [-2.969991788512712, -0.6604598778102144, 1.233811274632784], [-3.0567398040787133, -1.3216748454235834, -0.1856445039673828], [-3.2568593067482463, -2.6222799204493072, -0.06365857631325687], [-4.051038764929399, -0.7843632643516785, -0.8867620988473476], [-1.4799776082122547, -1.067353774287152, -1.170938927770444], [-1.4697169975946436, -2.0018705855472056, -2.242854995917617], [-1.3185419201178108, 0.32991570256916525, -1.335052605759722], [-0.4497836365109068, -1.541217683797412, -0.07791742489739739], [-2.2659422418234896, -1.4610539209006657, 2.0290187736421847], [-4.202897212097874, -0.5310943056708426, 1.7067156695461205], [-2.6212054717929303, 1.4189848017416486, 2.2776416618875777], [-2.7133779807411305, 1.3971249246106923, 0.11259346120335326], [-0.37104516182828423, 2.0617372872159114, 1.060259129110609], [-0.24552532002048982, 0.30567315029126085, 2.3240905646658514], [1.5770424436171642, 0.0, 0.0], [2.2927181468939155, 1.3915527243580577, 0.0], [3.782335574419715, 1.3186147352454591, -0.4807937764724355], [4.419592497958774, 0.34925073248439087, 0.14976318527134957], [3.8473329977166424, 1.1104432852891981, -1.7795740852228679], [4.36801041105158, 2.4684340853352054, -0.2056055442087168], [2.2938967322202157, 1.8704189044736064, 1.2405689893126328], [1.6494649440008806, 2.23521489431434, -0.8029305726284865], [1.997422457333483, -0.6906780683055267, 1.0535722235492984], [1.927718322430896, -0.6529932317206271, -1.1102241252095293], [-0.35014935725347385, -1.283813661620946, 0.08241309473865085], [-0.4266843221927557, 0.4915335256355437, -1.158605816601228], [-0.8284992432342889, -1.441086266844449, 0.8066788532199901]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0148', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
