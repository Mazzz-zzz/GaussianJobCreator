import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0086'
logfile = 'conf/5009017845242299296281_0086.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863816, 0.7718203945763835, 1.1635336229088475], [-2.2709622836291916, 0.7431123812655651, 1.1797556627389005], [-2.9970239643019005, 1.465687704195481, 2.368157397369682], [-4.4759542789311215, 1.8399334351693104, 2.004439902161414], [-5.036363599701665, 0.851319602488171, 1.329839260224275], [-5.178548597123875, 2.0769019739058328, 3.1085717675857167], [-4.5636267180482815, 3.3871491523350987, 0.9460919446874061], [-5.8470664672130335, 3.4362300625135536, 0.3355147732094933], [-4.008587173203819, 4.4411246791171415, 1.711874854036643], [-3.5016327596326593, 2.9784194514598163, -0.14278430802458875], [-3.0150768858013435, 0.6515380066656246, 3.419739973558244], [-2.351095681931019, 2.58257357437241, 2.6776285729604377], [-2.6427825570536148, 1.3431159399560075, 0.04548369333676697], [-2.667789377892033, -0.523270637705674, 1.1465664052870888], [-0.37104516182828207, 2.0617372872159128, 1.060259129110613], [-0.24552532002048844, 0.305673150291259, 2.3240905646658536], [1.577042443617167, 0.0, 0.0], [2.2927181468939164, 1.3915527243580572, 0.0], [1.600521547008245, 2.4407219045638975, -0.9357086002340277], [1.376069276137163, 1.9145350254105735, -2.125714308102145], [0.4611709133627535, 2.85194677898325, -0.41882897761344995], [2.4048649076934927, 3.47808001118291, -1.068930340330684], [3.5400592233304655, 1.2295174124846093, -0.4315210556927493], [2.308846803952254, 1.8960947387583758, 1.2310220414904687], [1.9974224573334842, -0.690678068305528, 1.053572223549295], [1.9277183224308956, -0.6529932317206251, -1.110224125209536], [-0.3501493572534771, -1.2838136616209432, 0.08241309473864497], [-0.4266843221927563, 0.4915335256355492, -1.1586058166012274], [-3.941591047707632, 2.5921028883669375, -0.9129522147373009]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0086', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
