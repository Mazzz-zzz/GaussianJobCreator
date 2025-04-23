import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0080'
logfile = 'conf/5009017845242299296281_0080.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863816, 0.7718203945763835, 1.1635336229088475], [-2.2709622836291916, 0.7431123812655651, 1.1797556627389005], [-2.9970239643019005, 1.465687704195481, 2.368157397369682], [-4.4759542789311215, 1.8399334351693104, 2.004439902161414], [-5.036363599701665, 0.851319602488171, 1.329839260224275], [-5.178548597123875, 2.0769019739058328, 3.1085717675857167], [-4.5636267180482815, 3.3871491523350987, 0.9460919446874061], [-5.8470664672130335, 3.4362300625135536, 0.3355147732094933], [-4.008587173203819, 4.4411246791171415, 1.711874854036643], [-3.5016327596326593, 2.9784194514598163, -0.14278430802458875], [-3.0150768858013435, 0.6515380066656246, 3.419739973558244], [-2.351095681931019, 2.58257357437241, 2.6776285729604377], [-2.6427825570536148, 1.3431159399560075, 0.04548369333676697], [-2.667789377892033, -0.523270637705674, 1.1465664052870888], [-0.37104516182828207, 2.0617372872159128, 1.060259129110613], [-0.24552532002048844, 0.305673150291259, 2.3240905646658536], [1.577042443617167, 0.0, 0.0], [2.2927181468939164, 1.3915527243580572, 0.0], [3.7823355744197142, 1.3186147352454598, -0.48079377647244853], [4.419592497958777, 0.3492507324843934, 0.14976318527135185], [3.8473329977166384, 1.1104432852892008, -1.779574085222876], [4.368010411051578, 2.4684340853352067, -0.20560554420872096], [2.293896732220219, 1.870418904473606, 1.2405689893126344], [1.649464944000881, 2.23521489431434, -0.8029305726284844], [1.9974224573334838, -0.6906780683055276, 1.0535722235492953], [1.9277183224308962, -0.6529932317206255, -1.1102241252095357], [-0.3501493572534771, -1.2838136616209432, 0.08241309473864497], [-0.4266843221927563, 0.4915335256355492, -1.1586058166012274], [-2.9999597154456636, 3.753921713140217, -0.4306780133408718]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0080', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
