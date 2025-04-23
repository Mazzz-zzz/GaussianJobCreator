import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0075'
logfile = 'conf/5009017845242299296281_0075.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863814, 0.6217394783082119, -1.2501828803165036], [-2.270962283629191, 0.6501421835576511, -1.2334320314121803], [-3.0203184489305857, -0.7271986321270443, -1.1788705129599806], [-3.131423238500184, -1.3851003159563176, -2.5981683179154134], [-4.113036047167058, -0.8202890821537163, -3.2792898247763516], [-1.9898347010198014, -1.2421684130768764, -3.2654359724481914], [-3.4853445173942257, -3.2247533095627285, -2.4882596168223086], [-4.470675829029669, -3.428858836708267, -1.48332471947352], [-3.5964475982228166, -3.7021802001725277, -3.8168214110949963], [-2.089155928253522, -3.6920001055840865, -1.9289024098694783], [-2.3454389866601137, -1.5526059487118047, -0.38374031881573706], [-4.247573342377323, -0.5532054861758633, -0.7055824261016228], [-2.596140664371264, 1.3384499968011874, -0.1355618835040283], [-2.689777026166647, 1.3195442638618216, -2.3005750802147893], [-0.3710451618282818, -0.11265730320381116, -2.315646431213898], [-0.24552532002048805, 1.859884894550714, -1.426765995739986], [1.5770424436171642, 0.0, 0.0], [2.2927181468939173, 1.3915527243580552, 0.0], [1.6005215470082501, 2.4407219045638957, -0.9357086002340295], [1.3760692761371698, 1.9145350254105753, -2.1257143081021446], [0.4611709133627533, 2.851946778983248, -0.41882897761345284], [2.404864907693497, 3.4780800111829073, -1.068930340330688], [3.540059223330467, 1.22951741248461, -0.4315210556927397], [2.3088468039522465, 1.8960947387583722, 1.2310220414904738], [1.997422457333482, -0.6906780683055239, 1.0535722235493052], [1.9277183224308936, -0.6529932317206333, -1.1102241252095246], [-0.35014935725347607, 0.7132786644586415, 1.0706086973199311], [-0.4266843221927586, -1.249148832966849, 0.15362238828850439], [-1.8506310768882523, -4.553881898325778, -2.297993262258889]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0075', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
