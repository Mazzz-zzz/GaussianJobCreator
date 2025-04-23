import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0137'
logfile = 'conf/5009017845242299296281_0137.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863861, 0.7718203945763877, 1.1635336229088449], [-0.3976197158559593, 2.318865504557597, 1.2514273698287335], [-1.1233882121466403, 3.1329845708321886, 2.379279124290352], [-2.5854106473305416, 3.5221161086622748, 1.9661784049894473], [-3.1611759047412464, 2.508175099936881, 1.3441495131679193], [-3.2994930518117376, 3.8458634914149767, 3.0405033057174085], [-2.6135892170537467, 5.00188486476642, 0.8124220495220797], [-2.3985135550467778, 6.174523331947342, 1.5877235990560843], [-1.846409481935056, 4.651091832232935, -0.3250032332974883], [-4.132502463634926, 4.940312367995723, 0.4005312790371705], [-1.1808982350936317, 2.3869269963298305, 3.4788618839288854], [-0.4526295328666409, 4.249846884503502, 2.630538289218348], [0.9186702629217367, 2.411095160260261, 1.4604772899922647], [-0.686339318432844, 2.8685142811653983, 0.07804489435392425], [-0.25604457595343416, 0.2684346684142553, 2.3183978473060134], [-2.007602477187449, 0.5917020341966386, 1.0837604470856856], [1.5770424436171646, 0.0, 0.0], [2.2927181468939146, 1.3915527243580523, 0.0], [2.341079856722322, 2.059852692894931, 1.416502376706472], [1.145391641007035, 2.0331725558915963, 1.9759511228307882], [3.202217992284638, 1.4422259961598234, 2.1984030628363325], [2.7215555591495044, 3.3152059779287946, 1.274535884539401], [1.6292323391939838, 2.212255867310577, -0.8090479336198905], [3.5455868300943827, 1.2600392214310707, -0.42809146886197624], [1.9974224573334807, -0.6906780683055224, 1.0535722235493095], [1.9277183224308982, -0.6529932317206333, -1.1102241252095253], [-0.35014935725347546, -1.2838136616209448, 0.08241309473865537], [-0.4266843221927533, 0.4915335256355349, -1.1586058166012312], [-4.241031985260075, 5.204613059938133, -0.5237678467173548]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0137', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
