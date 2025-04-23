import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0401'
logfile = 'conf/5009017845242299296281_0401.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586383, -1.393559872884599, 0.08664925740765088], [-2.2709622836291934, -1.3932545648232142, 0.05367636867327484], [-3.0203184489305848, -0.6573324959322145, 1.2192077454993062], [-4.480058160568046, -0.25598657285096993, 0.8097366435926605], [-5.037241623014006, -1.2341925500908644, 0.11752164928389947], [-5.212262690872639, -0.006028967204993568, 1.891539657886009], [-4.5062276883420225, 1.2924765630661823, -0.25009372190774054], [-4.31336332463281, 2.4163730003626984, 0.5996657400891442], [-3.71978584464473, 1.0199205179575954, -1.3957722110157407], [-6.0185340890925065, 1.243439481287459, -0.6872062299291296], [-3.086116569243988, -1.4709098648754981, 2.2693287026335534], [-2.3639288743924607, 0.44751541508172954, 1.5491980808953365], [-2.5961406643712643, -0.7866250332999412, -1.0913507571730168], [-2.6897770261666523, -2.6521285947103377, 0.007528686185012741], [-0.3710451618282815, -1.9490799840121078, 1.2553873021032877], [-0.24552532002049296, -2.165558044841983, -0.8973245689258771], [1.5770424436171662, 0.0, 0.0], [2.2927181468939177, 1.3915527243580554, 0.0], [1.6005215470082539, 2.440721904563898, -0.9357086002340279], [1.3760692761371682, 1.9145350254105777, -2.125714308102144], [0.4611709133627544, 2.8519467789832498, -0.4188289776134557], [2.404864907693502, 3.4780800111829073, -1.0689303403306858], [3.540059223330468, 1.2295174124846036, -0.43152105569274646], [2.308846803952251, 1.8960947387583706, 1.2310220414904662], [1.9974224573334842, -0.6906780683055239, 1.0535722235493021], [1.9277183224308954, -0.6529932317206275, -1.1102241252095344], [-0.3501493572534737, 0.5705349971623088, -1.1530217920585797], [-0.42668432219275576, 0.7576153073313022, 1.004983428312725], [-6.55502152370259, 1.7986803072302038, -0.1042386112675889]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0401', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
