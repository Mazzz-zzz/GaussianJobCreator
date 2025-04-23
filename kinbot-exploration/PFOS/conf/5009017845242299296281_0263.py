import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0263'
logfile = 'conf/5009017845242299296281_0263.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863819, -1.3935598728845997, 0.08664925740765092], [-2.270962283629191, -1.393254564823218, 0.053676368673273744], [-3.0203184489305834, -0.6573324959322201, 1.2192077454993062], [-4.480058160568046, -0.25598657285097753, 0.8097366435926645], [-5.202953426281216, -0.026515050940204468, 1.8919876304337158], [-4.464154537777423, 0.8343514247132503, 0.048283975068355564], [-5.325329553432012, -1.6109188884087284, -0.1758281811541025], [-6.725522199345814, -1.3624477614697956, -0.16343216499558497], [-4.563283221999027, -1.7983674024816865, -1.354600939139878], [-5.018129818548462, -2.818409682350801, 0.7875704566078667], [-3.0861165692439863, -1.4709098648755037, 2.269328702633553], [-2.3639288743924602, 0.4475154150817235, 1.549198080895336], [-2.5961406643712626, -0.7866250332999434, -1.0913507571730179], [-2.689777026166649, -2.6521285947103412, 0.00752868618500804], [-0.3710451618282798, -1.9490799840121118, 1.2553873021032853], [-0.2455253200204888, -2.165558044841983, -0.8973245689258781], [1.577042443617165, 0.0, 0.0], [2.292718146893917, 1.3915527243580557, 0.0], [1.600521547008257, 2.4407219045638993, -0.9357086002340261], [1.3760692761371667, 1.9145350254105789, -2.125714308102142], [0.46117091336275884, 2.8519467789832524, -0.4188289776134527], [2.404864907693507, 3.47808001118291, -1.0689303403306827], [3.5400592233304655, 1.229517412484604, -0.4315210556927458], [2.308846803952252, 1.8960947387583729, 1.2310220414904707], [1.997422457333484, -0.6906780683055261, 1.053572223549299], [1.9277183224308958, -0.6529932317206251, -1.1102241252095344], [-0.35014935725347524, 0.5705349971623078, -1.1530217920585801], [-0.4266843221927563, 0.7576153073313013, 1.0049834283127257], [-4.233362975845368, -2.6311306089335793, 1.3214524079321328]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0263', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
