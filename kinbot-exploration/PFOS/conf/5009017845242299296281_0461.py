import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0461'
logfile = 'conf/5009017845242299296281_0461.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863822, -1.3935598728845975, 0.08664925740765327], [-2.2709622836291916, -1.3932545648232157, 0.05367636867327366], [-2.9970239643018997, -2.7837283183799237, 0.08524408716293655], [-2.2704853996782592, -3.7973876868622507, 1.0361497385314784], [-3.0948301649158365, -4.7777601516461345, 1.3614610992245502], [-1.1948366276544546, -4.302842928203061, 0.43916166595565304], [-1.6927847351456848, -2.9753583691892023, 2.6211396672690737], [-2.7310271787076, -2.1200423017899976, 3.0825701181989356], [-1.0860228383589428, -3.983176550298476, 3.4094919274176454], [-0.5420771395216657, -2.0894271172534538, 2.011598466015423], [-4.239914773111563, -2.6098357638694947, 0.525670436701998], [-3.023054582933851, -3.305441500441076, -1.1344147028460818], [-2.6427825570536148, -0.710948003865581, 1.1404306775613293], [-2.667789377892029, -0.7313203152515887, -1.0264488679511448], [-0.37104516182828284, -1.9490799840121091, 1.2553873021032886], [-0.24552532002048746, -2.165558044841982, -0.8973245689258733], [1.5770424436171666, 0.0, 0.0], [2.292718146893919, 1.3915527243580557, 0.0], [2.341079856722325, 2.059852692894931, 1.4165023767064713], [1.1453916410070388, 2.0331725558916003, 1.975951122830787], [3.202217992284645, 1.4422259961598267, 2.198403062836327], [2.721555559149507, 3.3152059779287963, 1.2745358845394026], [1.629232339193983, 2.212255867310575, -0.8090479336198885], [3.5455868300943862, 1.2600392214310676, -0.4280914688619787], [1.9974224573334822, -0.6906780683055234, 1.0535722235493032], [1.9277183224308971, -0.6529932317206284, -1.1102241252095322], [-0.35014935725347157, 0.5705349971623073, -1.1530217920585843], [-0.42668432219275654, 0.7576153073313056, 1.004983428312722], [-0.8665600427703144, -1.1966315655545428, 1.8283245315294476]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0461', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
