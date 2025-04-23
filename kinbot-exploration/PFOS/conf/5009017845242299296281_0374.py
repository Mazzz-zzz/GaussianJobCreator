import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0374'
logfile = 'conf/5009017845242299296281_0374.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863814, 0.6217394783082185, -1.2501828803165003], [-2.270962283629191, 0.6501421835576563, -1.233432031412173], [-3.0203184489305848, -0.7271986321270386, -1.1788705129599821], [-2.2823447039971083, -1.7540796922695472, -0.2512417083665385], [-1.8454205729265496, -1.1434263302929446, 0.8362791963851384], [-3.1051738563474753, -2.7400657292245305, 0.0947117011053095], [-0.8077838749090355, -2.535088953994005, -1.1099373619800337], [-0.10865718900748936, -1.513764224017929, -1.8103291671931032], [-0.20471546093758583, -3.4117735479663707, -0.17554484145694824], [-1.592498989935659, -3.4164502728743646, -2.1529939489585317], [-4.245418844080691, -0.5314114745498586, -0.6995084752088649], [-3.094058093800262, -1.2477925357801305, -2.397055332165478], [-2.5961406643712626, 1.3384499968011894, -0.1355618835040201], [-2.6897770261666483, 1.3195442638618344, -2.3005750802147786], [-0.37104516182828196, -0.11265730320379871, -2.3156464312138967], [-0.2455253200204859, 1.8598848945507214, -1.4267659957399759], [1.5770424436171657, 0.0, 0.0], [2.2927181468939195, 1.3915527243580565, 0.0], [1.6005215470082523, 2.440721904563895, -0.9357086002340326], [1.3760692761371682, 1.914535025410577, -2.1257143081021455], [0.4611709133627604, 2.851946778983251, -0.4188289776134483], [2.4048649076935065, 3.478080011182904, -1.068930340330689], [3.5400592233304664, 1.2295174124846031, -0.43152105569274846], [2.3088468039522554, 1.8960947387583715, 1.2310220414904691], [1.997422457333483, -0.6906780683055278, 1.0535722235493001], [1.9277183224308947, -0.6529932317206282, -1.1102241252095306], [-0.3501493572534734, 0.7132786644586355, 1.0706086973199367], [-0.4266843221927593, -1.24914883296685, 0.15362238828850341], [-2.483578085755113, -3.607736409547624, -1.828417231361286]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0374', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
