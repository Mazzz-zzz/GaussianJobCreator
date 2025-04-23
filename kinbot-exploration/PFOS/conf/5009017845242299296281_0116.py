import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0116'
logfile = 'conf/5009017845242299296281_0116.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863845, -1.3935598728846021, 0.08664925740764984], [-0.39761971585595557, -2.243200645541641, 1.3824827499919123], [1.0879414097563276, -2.6761952078648314, 1.6417185410974557], [2.1006320219594956, -1.5474952427548012, 1.2416194637201219], [2.2762575910635965, -1.545484162421124, -0.06821012752820123], [1.6464483143989894, -0.3591354871283355, 1.6294574344592057], [3.7756073465590063, -1.7944622353315338, 2.0510194252196037], [3.6782290335556205, -1.3837659839366985, 3.40906105236394], [4.24869656716946, -3.0637074205497346, 1.6380302960762234], [4.578856613177446, -0.7022934487355706, 1.2494580819205137], [1.2395907752722297, -2.9383364670191585, 2.9368154419869623], [1.371737329774654, -3.7590542260525566, 0.9294643561329367], [-0.7843720924895153, -1.4713296274432253, 2.4020244202192758], [-1.158198376361052, -3.3308565421497485, 1.3537778366861704], [-0.2560445759534292, -2.142008766053285, -0.9267276815498303], [-2.0076024771874486, -1.234415095891312, -0.029451230457635254], [1.5770424436171642, 0.0, 0.0], [2.2927181468939124, 1.3915527243580568, 0.0], [1.6005215470082428, 2.4407219045638975, -0.9357086002340265], [1.376069276137152, 1.9145350254105755, -2.125714308102144], [0.4611709133627475, 2.8519467789832524, -0.41882897761344284], [2.4048649076934936, 3.47808001118291, -1.0689303403306885], [3.5400592233304575, 1.2295174124846064, -0.43152105569275384], [2.308846803952254, 1.8960947387583766, 1.2310220414904651], [1.9974224573334824, -0.6906780683055287, 1.0535722235492933], [1.9277183224308927, -0.6529932317206194, -1.1102241252095364], [-0.35014935725347973, 0.5705349971623124, -1.1530217920585764], [-0.4266843221927598, 0.7576153073313006, 1.0049834283127321], [5.484510059973413, -1.0019713820084917, 1.0883762477123806]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0116', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
