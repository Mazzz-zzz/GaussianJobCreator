import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0124'
logfile = 'conf/5009017845242299296281_0124.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586383, -1.393559872884599, 0.08664925740765088], [-0.397619715855956, -2.243200645541635, 1.38248274999192], [1.0879414097563267, -2.6761952078648257, 1.6417185410974653], [1.3412747433709995, -3.0081801923381053, 3.1534041404090263], [1.4942413294616161, -1.8890114354210112, 3.839381151800391], [0.3200544491791535, -3.6981149368927495, 3.6534033399484347], [2.8849342358987173, -4.048842112015683, 3.3896274906033725], [2.5770387278690783, -5.389519592485756, 3.0287739783500136], [3.9655135290297987, -3.3100663187062995, 2.849335449088397], [2.9826164228846865, -3.9579511190762653, 4.958939329454128], [1.3507067731889262, -3.7633797193128684, 0.9220132164307552], [1.9101736356829004, -1.6963410126419183, 1.2892277976148896], [-0.7843720924895166, -1.471329627443213, 2.4020244202192753], [-1.1581983763610533, -3.3308565421497396, 1.3537778366861806], [-0.2560445759534264, -2.1420087660532863, -0.926727681549823], [-2.007602477187446, -1.2344150958913094, -0.029451230457633544], [1.577042443617165, 0.0, 0.0], [2.2927181468939164, 1.3915527243580577, 0.0], [1.6005215470082519, 2.4407219045638993, -0.9357086002340291], [1.3760692761371667, 1.914535025410577, -2.125714308102147], [0.4611709133627542, 2.8519467789832547, -0.4188289776134507], [2.404864907693501, 3.478080011182911, -1.0689303403306871], [3.540059223330464, 1.2295174124846064, -0.43152105569274873], [2.308846803952253, 1.8960947387583742, 1.231022041490468], [1.9974224573334836, -0.6906780683055249, 1.0535722235493], [1.927718322430893, -0.652993231720625, -1.110224125209532], [-0.3501493572534752, 0.5705349971623124, -1.1530217920585777], [-0.4266843221927579, 0.7576153073313051, 1.004983428312727], [2.565677221337377, -4.72978031493891, 5.366865525380696]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0124', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
