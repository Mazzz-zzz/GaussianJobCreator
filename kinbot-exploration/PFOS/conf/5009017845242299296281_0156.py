import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0156'
logfile = 'conf/5009017845242299296281_0156.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863814, -1.393559872884599, 0.08664925740765088], [-2.270962283629192, -1.3932545648232202, 0.053676368673273835], [-3.0203184489305857, -0.6573324959322203, 1.2192077454993038], [-3.131423238500184, -1.5575296086444828, 2.4986162193657244], [-4.113036047167058, -2.4298037535512926, 2.3500360959803044], [-1.9898347010198019, -2.2068663000332487, 2.7084673877272656], [-3.4853445173942257, -0.5425193845977042, 4.036848095430419], [-3.9613100631246523, -1.420500303166332, 5.049289359759823], [-2.39398443908955, 0.34415144242193874, 4.204492945166777], [-4.698972842349626, 0.2945479739521012, 3.48280192774634], [-2.3454389866601146, 0.44397410980512847, 1.536466353059133], [-4.247573342377323, -0.3344495623799341, 0.8318812175920285], [-2.596140664371265, -0.7866250332999446, -1.0913507571730185], [-2.6897770261666474, -2.6521285947103417, 0.0075286861850103965], [-0.37104516182828184, -1.949079984012111, 1.255387302103287], [-0.24552532002048807, -2.165558044841982, -0.8973245689258756], [1.5770424436171642, 0.0, 0.0], [2.292718146893914, 1.391552724358057, 0.0], [1.6005215470082497, 2.4407219045639006, -0.9357086002340282], [1.3760692761371678, 1.9145350254105806, -2.1257143081021406], [0.46117091336275107, 2.851946778983253, -0.41882897761345117], [2.404864907693502, 3.4780800111829127, -1.0689303403306858], [3.5400592233304646, 1.229517412484608, -0.4315210556927457], [2.30884680395225, 1.896094738758375, 1.2310220414904705], [1.997422457333483, -0.6906780683055248, 1.0535722235493], [1.9277183224308945, -0.6529932317206222, -1.110224125209533], [-0.3501493572534773, 0.570534997162307, -1.1530217920585837], [-0.4266843221927594, 0.7576153073313044, 1.0049834283127241], [-4.6564733123172175, 0.3497911583921597, 2.5178644338038865]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0156', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
