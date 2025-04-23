import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0490'
logfile = 'conf/5009017845242299296281_0490.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863822, 0.6217394783082169, -1.2501828803164996], [-2.2709622836291916, 0.6501421835576559, -1.2334320314121736], [-2.9970239643018997, 1.3180406141844494, -2.4534014845326153], [-2.2704853996782592, 1.0013618477382709, -3.806709074106673], [-3.0948301649158365, 1.209820177630333, -4.818392214126816], [-1.1948366276544546, 1.7710963050156465, -3.945952117295896], [-1.6927847351456848, -0.7822943541274964, -3.8873057667150266], [-2.7310271787076, -1.6095628804120714, -3.3772955495472448], [-1.0860228383589428, -0.9611183479923998, -5.15427804402577], [-0.5420771395216657, -0.6973818151564307, -2.8152961959052907], [-4.239914773111563, 0.8496739297323636, -2.523019289567145], [-3.0230545829338507, 2.6351527013118248, -2.2953889586822775], [-2.6427825570536148, -0.632167936090418, -1.1859143708980966], [-2.667789377892029, 1.2545909529572645, -0.12011753733594527], [-0.37104516182828284, -0.1126573032037983, -2.315646431213899], [-0.24552532002048746, 1.8598848945507207, -1.4267659957399768], [1.5770424436171666, 0.0, 0.0], [2.29271814689392, 1.3915527243580577, 0.0], [1.600521547008252, 2.4407219045639, -0.9357086002340262], [1.376069276137165, 1.914535025410578, -2.125714308102147], [0.4611709133627617, 2.851946778983253, -0.4188289776134476], [2.4048649076935034, 3.478080011182911, -1.0689303403306845], [3.540059223330466, 1.229517412484606, -0.43152105569274923], [2.308846803952257, 1.8960947387583726, 1.2310220414904716], [1.997422457333483, -0.6906780683055281, 1.053572223549299], [1.9277183224308936, -0.652993231720624, -1.1102241252095326], [-0.3501493572534734, 0.7132786644586327, 1.0706086973199345], [-0.42668432219275654, -1.2491488329668474, 0.15362238828850197], [-0.8665600427703144, -0.9850597078895069, -1.9504756005053048]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0490', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
