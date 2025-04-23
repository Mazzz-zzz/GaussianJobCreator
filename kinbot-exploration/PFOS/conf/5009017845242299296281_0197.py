import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0197'
logfile = 'conf/5009017845242299296281_0197.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863809, -1.393559872884599, 0.08664925740765585], [-0.3466020415139007, -2.420439425248669, -1.059551387211253], [-0.7363023803695176, -2.0273092602548024, -2.5275982746061025], [-2.245266873553548, -2.329948061955714, -2.8292952633496165], [-2.4038439949868793, -3.6181636043108996, -3.0780038399396776], [-3.0033804050952075, -1.9822263526669373, -1.7933157622320852], [-2.8645095706331407, -1.3711226263309546, -4.31887761999445], [-4.070659490453873, -1.976453114053919, -4.7674401901666705], [-2.7357825244594305, 0.003827407378341697, -4.005351665671368], [-1.7096109057685829, -1.7411654311169051, -5.323775637869506], [-0.5266683954723095, -0.7244611301521561, -2.693826716126075], [0.004125934750529869, -2.71341613896585, -3.388674962621095], [0.9809271994195111, -2.562327183726283, -1.0109564779824374], [-0.9145558447621667, -3.5839764719224845, -0.7664486088494463], [-2.0119541879597156, -1.210216503485702, 0.0009668077396042128], [-0.40994706586377383, -1.949605403411485, 1.259403781369357], [1.5770424436171668, 0.0, 0.0], [2.2927181468939093, 1.3915527243580605, 0.0], [3.7823355744197142, 1.318614735245461, -0.48079377647244564], [4.419592497958774, 0.3492507324844044, 0.14976318527134708], [3.8473329977166353, 1.1104432852891999, -1.7795740852228807], [4.368010411051576, 2.468434085335212, -0.20560554420872262], [2.2938967322202135, 1.8704189044736106, 1.240568989312631], [1.6494649440008806, 2.235214894314338, -0.8029305726284872], [1.9974224573334851, -0.6906780683055231, 1.0535722235493004], [1.927718322430895, -0.6529932317206257, -1.110224125209531], [-0.3501493572534788, 0.5705349971623067, -1.1530217920585792], [-0.4266843221927566, 0.7576153073313048, 1.0049834283127286], [-1.514435757135947, -0.9908431528998533, -5.902468584831365]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0197', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
