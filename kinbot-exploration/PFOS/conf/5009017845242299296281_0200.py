import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0200'
logfile = 'conf/5009017845242299296281_0200.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863756, 0.621739478308203, -1.2501828803165123], [-2.2709622836291863, 0.6501421835576423, -1.2334320314121954], [-3.0203184489305785, -0.727198632127055, -1.1788705129599912], [-3.131423238500171, -1.3851003159563369, -2.598168317915419], [-4.113036047167041, -0.8202890821537401, -3.279289824776366], [-1.9898347010197868, -1.242168413076899, -3.2654359724481945], [-3.485344517394217, -3.2247533095627445, -2.4882596168223006], [-2.275694716609433, -3.8908785884863324, -2.1484722126262312], [-4.69823536319325, -3.368017453425294, -1.7715114372583425], [-3.796615807083963, -3.487819958394286, -4.009590860272073], [-2.3454389866601097, -1.5526059487118096, -0.3837403188157383], [-4.247573342377318, -0.553205486175873, -0.705582426101639], [-2.5961406643712612, 1.338449996801184, -0.13556188350404783], [-2.689777026166637, 1.3195442638618058, -2.3005750802148097], [-0.3710451618282712, -0.1126573032038257, -2.3156464312139007], [-0.24552532002048144, 1.859884894550705, -1.4267659957400014], [1.5770424436171677, 0.0, 0.0], [2.292718146893916, 1.3915527243580565, 0.0], [1.600521547008243, 2.440721904563891, -0.9357086002340365], [1.3760692761371716, 1.9145350254105655, -2.1257143081021503], [0.46117091336273863, 2.851946778983237, -0.4188289776134648], [2.404864907693477, 3.478080011182911, -1.0689303403306925], [3.540059223330462, 1.2295174124846244, -0.4315210556927345], [2.30884680395223, 1.8960947387583813, 1.2310220414904722], [1.9974224573334813, -0.6906780683055174, 1.0535722235493121], [1.927718322430899, -0.6529932317206386, -1.1102241252095175], [-0.3501493572534797, 0.7132786644586477, 1.0706086973199238], [-0.4266843221927574, -1.2491488329668485, 0.15362238828850885], [-4.493119589630578, -4.153163286648216, -4.0999408638262285]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0200', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
