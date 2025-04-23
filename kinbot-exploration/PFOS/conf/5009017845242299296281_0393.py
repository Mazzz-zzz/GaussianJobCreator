import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0393'
logfile = 'conf/5009017845242299296281_0393.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863817, -1.393559872884598, 0.08664925740765204], [-2.2709622836291916, -1.3932545648232184, 0.05367636867327128], [-3.020318448930586, -0.6573324959322178, 1.2192077454993018], [-3.131423238500186, -1.557529608644479, 2.4986162193657226], [-3.3784464002921473, -0.8053658014461091, 3.556789696521597], [-4.1058597005417825, -2.45093059767987, 2.3529855399487376], [-1.5468226482012561, -2.509766594543162, 2.821027483166894], [-1.5707570250137253, -2.9706819695875013, 4.166172520129203], [-1.3349733932555896, -3.3474160961895847, 1.6990560141809483], [-0.534085487627692, -1.3080772354250902, 2.716870113841465], [-2.345438986660115, 0.44397410980513186, 1.5364663530591292], [-4.247573342377323, -0.3344495623799334, 0.8318812175920287], [-2.5961406643712617, -0.786625033299947, -1.0913507571730219], [-2.689777026166647, -2.6521285947103412, 0.007528686185010396], [-0.3710451618282826, -1.949079984012108, 1.255387302103288], [-0.24552532002048646, -2.1655580448419816, -0.8973245689258742], [1.5770424436171642, 0.0, 0.0], [2.2927181468939137, 1.3915527243580597, 0.0], [1.6005215470082503, 2.440721904563901, -0.9357086002340304], [1.3760692761371678, 1.9145350254105777, -2.1257143081021446], [0.4611709133627495, 2.851946778983252, -0.41882897761345605], [2.4048649076934967, 3.478080011182912, -1.068930340330688], [3.540059223330463, 1.2295174124846124, -0.4315210556927463], [2.308846803952248, 1.896094738758379, 1.2310220414904693], [1.9974224573334811, -0.6906780683055214, 1.0535722235493041], [1.9277183224308956, -0.6529932317206221, -1.1102241252095302], [-0.35014935725347546, 0.5705349971623076, -1.153021792058585], [-0.42668432219276026, 0.7576153073313056, 1.004983428312722], [-0.33954832241438326, -0.9553362554214117, 3.596467453010255]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0393', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
