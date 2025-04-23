import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0272'
logfile = 'conf/5009017845242299296281_0272.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863789, -1.3935598728846015, 0.08664925740765352], [-2.2709622836291903, -1.3932545648232202, 0.05367636867327507], [-3.0203184489305848, -0.65733249593222, 1.2192077454993033], [-2.2823447039971083, 0.6594581441991471, 1.644698427951089], [-1.2544029753924617, 0.36695518978720976, 2.42205276628759], [-1.8529952461325305, 1.317133480284267, 0.5713754402676374], [-3.414407898966413, 1.8184222874767133, 2.591745611732225], [-2.6069712798285622, 2.7693497180464384, 3.2745443636901546], [-4.476172266748106, 2.1600544728539055, 1.719220168623798], [-3.9728808110846567, 0.7844262822367575, 3.6403174487500727], [-4.24541884408069, -0.34008637241846923, 0.8099700744271578], [-3.0940580938002618, -1.452014544042184, 2.279147700720941], [-2.5961406643712612, -0.7866250332999513, -1.0913507571730185], [-2.689777026166647, -2.652128594710344, 0.007528686185015114], [-0.3710451618282796, -1.9490799840121105, 1.2553873021032917], [-0.2455253200204863, -2.1655580448419856, -0.8973245689258746], [1.5770424436171668, 0.0, 0.0], [2.2927181468939155, 1.3915527243580548, 0.0], [1.6005215470082548, 2.440721904563895, -0.9357086002340316], [1.376069276137169, 1.9145350254105744, -2.125714308102145], [0.4611709133627544, 2.8519467789832507, -0.41882897761345117], [2.404864907693506, 3.4780800111829118, -1.068930340330682], [3.5400592233304655, 1.2295174124846084, -0.4315210556927432], [2.3088468039522496, 1.8960947387583749, 1.231022041490473], [1.9974224573334838, -0.6906780683055231, 1.0535722235493057], [1.9277183224308962, -0.6529932317206254, -1.1102241252095295], [-0.3501493572534736, 0.5705349971623042, -1.1530217920585832], [-0.42668432219275926, 0.7576153073313027, 1.0049834283127217], [-3.4641170317638257, 0.8308976579202422, 4.461878289339052]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0272', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
