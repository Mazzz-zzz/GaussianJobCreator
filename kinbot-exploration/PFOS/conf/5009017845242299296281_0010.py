import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0010'
logfile = 'conf/5009017845242299296281_0010.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863814, -1.3935598728845988, 0.0866492574076521], [-2.270962283629191, -1.3932545648232195, 0.05367636867327381], [-3.0203184489305857, -0.6573324959322182, 1.2192077454993027], [-3.131423238500184, -1.5575296086444828, 2.4986162193657244], [-4.113036047167058, -2.4298037535512904, 2.350036095980306], [-1.9898347010198014, -2.2068663000332482, 2.7084673877272647], [-3.4853445173942257, -0.5425193845977008, 4.036848095430419], [-4.470675829029668, 0.429832529228626, 3.7111412183168793], [-3.5964475982228166, -1.4543742036303895, 5.11459280828466], [-2.089155928253522, 0.1755215644240368, 4.161817087145389], [-2.3454389866601137, 0.44397410980512836, 1.5364663530591325], [-4.247573342377323, -0.3344495623799334, 0.8318812175920287], [-2.596140664371264, -0.7866250332999447, -1.0913507571730179], [-2.689777026166647, -2.6521285947103412, 0.007528686185012751], [-0.3710451618282818, -1.9490799840121091, 1.2553873021032886], [-0.24552532002048805, -2.1655580448419816, -0.8973245689258754], [1.5770424436171642, 0.0, 0.0], [2.292718146893915, 1.3915527243580559, 0.0], [1.6005215470082517, 2.440721904563902, -0.9357086002340262], [1.3760692761371698, 1.9145350254105828, -2.125714308102141], [0.46117091336275395, 2.8519467789832555, -0.4188289776134474], [2.404864907693505, 3.4780800111829144, -1.0689303403306798], [3.5400592233304655, 1.2295174124846082, -0.43152105569274235], [2.30884680395225, 1.8960947387583738, 1.231022041490472], [1.9974224573334818, -0.690678068305523, 1.0535722235493], [1.9277183224308938, -0.6529932317206224, -1.1102241252095328], [-0.35014935725347673, 0.5705349971623067, -1.1530217920585843], [-0.42668432219275976, 0.7576153073313053, 1.0049834283127235], [-1.8506310768882523, 0.28682040632119976, 5.092774040913673]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0010', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
