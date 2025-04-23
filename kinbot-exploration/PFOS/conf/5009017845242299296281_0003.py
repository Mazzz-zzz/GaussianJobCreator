import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0003'
logfile = 'conf/5009017845242299296281_0003.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863822, 0.6217394783082169, -1.2501828803164996], [-2.2709622836291916, 0.6501421835576559, -1.2334320314121736], [-2.9970239643018997, 1.3180406141844494, -2.4534014845326153], [-2.2704853996782592, 1.0013618477382709, -3.806709074106673], [-3.0948301649158365, 1.209820177630333, -4.818392214126816], [-1.1948366276544546, 1.7710963050156465, -3.945952117295896], [-1.6927847351456848, -0.7822943541274964, -3.8873057667150266], [-2.7310271787076, -1.6095628804120714, -3.3772955495472448], [-1.0860228383589428, -0.9611183479923998, -5.15427804402577], [-0.5420771395216657, -0.6973818151564307, -2.8152961959052907], [-4.239914773111563, 0.8496739297323636, -2.523019289567145], [-3.0230545829338507, 2.6351527013118248, -2.2953889586822775], [-2.6427825570536148, -0.632167936090418, -1.1859143708980966], [-2.667789377892029, 1.2545909529572645, -0.12011753733594527], [-0.37104516182828284, -0.1126573032037983, -2.315646431213899], [-0.24552532002048746, 1.8598848945507207, -1.4267659957399768], [1.5770424436171666, 0.0, 0.0], [2.29271814689392, 1.3915527243580577, 0.0], [2.3410798567223297, 2.059852692894926, 1.416502376706475], [1.1453916410070395, 2.033172555891594, 1.975951122830798], [3.2022179922846488, 1.4422259961598265, 2.198403062836323], [2.7215555591495075, 3.315205977928796, 1.274535884539404], [1.6292323391939771, 2.212255867310578, -0.809047933619883], [3.5455868300943836, 1.2600392214310738, -0.4280914688619856], [1.9974224573334833, -0.6906780683055279, 1.0535722235492995], [1.9277183224308934, -0.6529932317206238, -1.110224125209533], [-0.3501493572534734, 0.7132786644586327, 1.0706086973199345], [-0.42668432219275654, -1.2491488329668474, 0.15362238828850197], [0.19841530887108264, -1.2648388342874792, -3.071471627998648]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0003', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
