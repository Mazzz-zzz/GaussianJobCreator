import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0041'
logfile = 'conf/5009017845242299296281_0041.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863838, -1.3935598728845957, 0.08664925740765067], [-2.270962283629194, -1.3932545648232146, 0.05367636867327237], [-2.9970239643019028, -2.7837283183799224, 0.08524408716293154], [-3.0516123514917357, -3.443374587786912, -1.3365407299748289], [-3.279316018624865, -4.740067389005451, -1.2214722185603748], [-4.013520036470852, -2.8879361664696774, -2.0680609321946637], [-1.440485817204604, -3.2181908870364535, -2.272076364008747], [-0.3660848203455984, -3.456054939018342, -1.3712355353752257], [-1.5924589367306476, -3.875980277254094, -3.5168952507174285], [-1.5450502463939484, -1.6689398312376449, -2.5355722388599236], [-2.3329497198238234, -3.597057518120825, 0.9016340449500758], [-4.241488402883606, -2.6316641353094203, 0.51954185001882], [-2.6427825570536183, -0.7109480038655804, 1.1404306775613284], [-2.6677893778920336, -0.7313203152515851, -1.0264488679511463], [-0.3710451618282875, -1.949079984012109, 1.2553873021032884], [-0.24552532002048894, -2.16555804484198, -0.8973245689258756], [1.5770424436171664, 0.0, 0.0], [2.292718146893916, 1.391552724358056, 0.0], [3.7823355744197147, 1.3186147352454594, -0.48079377647243915], [4.41959249795878, 0.3492507324843911, 0.14976318527135712], [3.8473329977166397, 1.1104432852891952, -1.779574085222881], [4.368010411051583, 2.468434085335202, -0.20560554420871308], [2.293896732220213, 1.8704189044736095, 1.2405689893126335], [1.6494649440008873, 2.2352148943143346, -0.8029305726284868], [1.9974224573334805, -0.6906780683055245, 1.053572223549303], [1.9277183224308947, -0.6529932317206246, -1.1102241252095315], [-0.35014935725347235, 0.5705349971623097, -1.1530217920585826], [-0.4266843221927606, 0.7576153073313081, 1.0049834283127235], [-1.1714428886352113, -1.4491909608674787, -3.40049441044276]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0041', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
