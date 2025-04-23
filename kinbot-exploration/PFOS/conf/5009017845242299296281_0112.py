import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0112'
logfile = 'conf/5009017845242299296281_0112.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863861, 0.7718203945763867, 1.1635336229088478], [-0.39761971585595707, 2.3188655045575937, 1.2514273698287424], [1.0879414097563254, 2.759867566386749, 1.4967937649483767], [1.4730420899455337, 2.6857139528931335, 3.015274540831185], [2.7879690257787098, 2.6561530893205583, 3.144097100354801], [0.9875700214111517, 3.7379293804663556, 3.6680069555333064], [0.7750378848301906, 1.150016923244951, 3.8375203032927736], [0.9638915411618125, 0.04927612753627226, 2.9571027858825403], [1.2095491862312482, 1.1716267414164387, 5.185144032811274], [-0.746936443769862, 1.5548071048149876, 3.8207369247317544], [1.8947536492841401, 1.9494330700686107, 0.8176461448894714], [1.259810370563855, 4.009816458287681, 1.086283558698188], [-0.7843720924895176, 2.81587898214209, 0.07319662459686953], [-1.1581983763610557, 2.8378342687254494, 2.2077174635201686], [-0.25604457595343133, 0.26843466841424685, 2.3183978473060116], [-2.007602477187449, 0.591702034196638, 1.0837604470856914], [1.577042443617165, 0.0, 0.0], [2.2927181468939164, 1.3915527243580545, 0.0], [2.341079856722327, 2.0598526928949292, 1.416502376706474], [1.1453916410070353, 2.0331725558915945, 1.9759511228307916], [3.2022179922846483, 1.4422259961598267, 2.198403062836325], [2.7215555591495075, 3.315205977928793, 1.2745358845394046], [1.6292323391939805, 2.212255867310575, -0.8090479336198871], [3.5455868300943836, 1.2600392214310687, -0.42809146886198146], [1.997422457333483, -0.6906780683055265, 1.0535722235493008], [1.927718322430891, -0.6529932317206273, -1.1102241252095326], [-0.350149357253475, -1.2838136616209441, 0.08241309473865074], [-0.4266843221927582, 0.49153352563554115, -1.15860581660123], [-0.8383239368037281, 2.517287628191143, 3.785495543735128]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0112', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
