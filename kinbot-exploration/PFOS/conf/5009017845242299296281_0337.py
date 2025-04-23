import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0337'
logfile = 'conf/5009017845242299296281_0337.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.69372834458638, -1.3935598728845973, 0.08664925740765574], [-0.3466020415139041, -2.4204394252486643, -1.059551387211261], [1.162453549046706, -2.799902192493109, -1.2603960465750481], [1.9172819807099777, -2.925959600137113, 0.1085191179741827], [2.227217092188635, -1.7242241351689576, 0.5626931396745398], [1.1566135117755332, -3.5540341854541517, 1.0005986190281884], [3.5062602700141103, -3.9101463550175772, -0.059372986524422224], [4.313987042727912, -3.6487677738874984, 1.0815431829541988], [3.1421067329679917, -5.208090461909852, -0.4930665056003186], [4.116322039547438, -3.1481171623606428, -1.2953660191510206], [1.2315520144708516, -3.970118415502823, -1.8887358797759217], [1.768541893208141, -1.8669832707228684, -1.983409881825592], [-1.0117971004258957, -3.5349698356523596, -0.7430628651027303], [-0.8215282967123818, -1.946695641113878, -2.205112865005799], [-2.011954187959715, -1.2102165034856955, 0.000966807739600983], [-0.40994706586377555, -1.9496054034114851, 1.2594037813693522], [1.5770424436171662, 0.0, 0.0], [2.2927181468939173, 1.391552724358055, 0.0], [2.341079856722328, 2.0598526928949292, 1.4165023767064713], [1.145391641007046, 2.0331725558915923, 1.9759511228307978], [3.2022179922846536, 1.442225996159827, 2.1984030628363262], [2.7215555591495066, 3.3152059779287955, 1.274535884539403], [1.6292323391939765, 2.212255867310575, -0.809047933619886], [3.545586830094382, 1.2600392214310745, -0.42809146886198857], [1.9974224573334842, -0.6906780683055272, 1.0535722235492957], [1.9277183224308938, -0.6529932317206232, -1.1102241252095353], [-0.3501493572534749, 0.5705349971623113, -1.1530217920585795], [-0.426684322192755, 0.7576153073313029, 1.0049834283127277], [3.4139916343943835, -2.727618423126553, -1.8109968939016334]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0337', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
