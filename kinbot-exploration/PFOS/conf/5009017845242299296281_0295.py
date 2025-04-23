import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0295'
logfile = 'conf/5009017845242299296281_0295.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863843, -1.3935598728845973, 0.08664925740764828], [-2.270962283629195, -1.393254564823212, 0.053676368673273515], [-2.9699917885127136, -0.7382819684025383, -1.1888806696804055], [-3.0567398040787124, 0.821610279220512, -1.0517817396960003], [-1.9299416078656424, 1.2912100611157198, -0.5455171715777618], [-3.2751648795012236, 1.3793631286707715, -2.2392590177564293], [-4.462079523721016, 1.345419937332617, 0.07618414878608065], [-4.239230824117428, 2.694401446526589, 0.467359752253039], [-5.660645840812834, 0.8789292575805179, -0.5164219606031193], [-4.136425722847151, 0.4051862494257172, 1.2970251883610098], [-2.2659422418234905, -1.0266548422793398, -2.2798191986199305], [-4.2028972120978745, -1.2125119740284827, -1.3132989952892684], [-2.621205471792931, -2.6819879407832747, 0.0900560549484987], [-2.71337798074113, -0.7960712600074682, 1.1536489463716026], [-0.37104516182828434, -1.9490799840121111, 1.2553873021032849], [-0.2455253200204931, -2.165558044841977, -0.8973245689258814], [1.5770424436171637, 0.0, 0.0], [2.2927181468939173, 1.3915527243580523, 0.0], [2.341079856722335, 2.0598526928949314, 1.416502376706472], [1.1453916410070453, 2.0331725558916, 1.9759511228307918], [3.2022179922846528, 1.442225996159828, 2.198403062836322], [2.721555559149511, 3.3152059779287963, 1.2745358845393975], [1.6292323391939747, 2.212255867310576, -0.8090479336198853], [3.5455868300943827, 1.2600392214310703, -0.4280914688619881], [1.997422457333484, -0.6906780683055276, 1.0535722235492933], [1.9277183224308925, -0.6529932317206237, -1.1102241252095348], [-0.3501493572534758, 0.5705349971623134, -1.1530217920585786], [-0.4266843221927538, 0.7576153073313038, 1.004983428312727], [-3.6343716860935866, -0.3672797694330289, 1.0017280827470718]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0295', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
