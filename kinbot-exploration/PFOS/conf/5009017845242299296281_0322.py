import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0322'
logfile = 'conf/5009017845242299296281_0322.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863829, 0.7718203945763835, 1.1635336229088495], [-2.2709622836291934, 0.7431123812655637, 1.179755662738898], [-3.0203184489305848, 1.384531128059259, -0.040337232539325606], [-4.480058160568045, 0.8292457901518746, -0.18317744667967567], [-5.037241623014006, 0.7188730088199354, 1.0100812768982363], [-5.212262690872638, 1.6411358795975055, -0.9405485901848986], [-4.5062276883420225, -0.8628257980321968, -0.9942706764574425], [-4.31336332463281, -0.6888607354849567, -2.392473273377493], [-3.71978584464473, -1.7187344516148015, -0.18539097288438794], [-6.018534089092507, -1.2168577934012865, -0.7332470638989184], [-3.086116569243988, 2.7007512384555885, 0.13918095834253655], [-2.363928874392461, 1.1178871860085953, -1.162158758493585], [-2.5961406643712643, -0.551824963501243, 1.2269126406770379], [-2.6897770261666523, 1.3325843308485132, 2.2930463940297683], [-0.3710451618282815, 2.061737287215912, 1.0602591291106067], [-0.24552532002049293, 0.3056731502912636, 2.324090564665856], [1.5770424436171662, 0.0, 0.0], [2.2927181468939173, 1.3915527243580543, 0.0], [1.6005215470082557, 2.4407219045638975, -0.9357086002340284], [1.3760692761371687, 1.9145350254105749, -2.1257143081021423], [0.4611709133627593, 2.8519467789832507, -0.41882897761345206], [2.404864907693505, 3.4780800111829118, -1.0689303403306827], [3.54005922333047, 1.2295174124846024, -0.43152105569274196], [2.308846803952252, 1.8960947387583715, 1.2310220414904722], [1.9974224573334844, -0.6906780683055244, 1.053572223549302], [1.9277183224308976, -0.6529932317206242, -1.1102241252095284], [-0.3501493572534737, -1.2838136616209428, 0.08241309473865066], [-0.42668432219275576, 0.4915335256355448, -1.1586058166012259], [-6.352055945044676, -0.7245415440417452, 0.02987171382495119]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0322', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
