import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0225'
logfile = 'conf/5009017845242299296281_0225.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863809, 0.6217394783082097, -1.2501828803165052], [-0.3466020415139007, 2.127818130564319, -1.5663863369811313], [-0.7363023803695176, 3.2026189464980006, -0.49190218340505043], [-2.2452668735535477, 3.6152156038456074, -0.6031465794771648], [-2.4038439949868793, 4.474711320489259, -1.594419676411675], [-3.003380405095207, 2.5441701834335073, -0.8200004963445009], [-2.8645095706331407, 4.4258190479167485, 0.9720117838909723], [-4.070659490453872, 5.116950872734213, 0.6720614889237702], [-2.7357825244594305, 3.4668225898725518, 2.0059904648559552], [-1.7096109057685829, 5.481107662002151, 1.1539943233962209], [-0.5266683954723095, 2.695152934634472, 0.7195116152968809], [0.004125934750529869, 4.291386672281075, -0.6555498260725717], [0.980927199419511, 2.1566775839163737, -1.713562194923181], [-0.9145558447621667, 2.4557522019201006, -2.7205903668258777], [-2.0119541879597156, 0.604270971679777, -1.0485616399675994], [-0.40994706586377383, -0.11587296658231004, -2.3181096973944326], [1.5770424436171668, 0.0, 0.0], [2.2927181468939146, 1.3915527243580554, 0.0], [1.600521547008242, 2.440721904563893, -0.9357086002340342], [1.3760692761371633, 1.914535025410568, -2.1257143081021477], [0.4611709133627453, 2.8519467789832467, -0.41882897761345417], [2.4048649076934883, 3.4780800111829064, -1.0689303403306913], [3.540059223330465, 1.2295174124846087, -0.4315210556927449], [2.3088468039522434, 1.8960947387583773, 1.2310220414904698], [1.9974224573334818, -0.6906780683055267, 1.0535722235493041], [1.9277183224308951, -0.6529932317206316, -1.1102241252095257], [-0.3501493572534771, 0.7132786644586372, 1.070608697319929], [-0.426684322192757, -1.2491488329668508, 0.15362238828850236], [-1.514435757135947, 5.607109315953477, 2.093138950838532]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0225', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
