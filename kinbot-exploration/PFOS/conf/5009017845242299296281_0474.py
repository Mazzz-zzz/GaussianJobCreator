import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0474'
logfile = 'conf/5009017845242299296281_0474.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863824, -1.393559872884598, 0.0866492574076533], [-2.2709622836291925, -1.3932545648232186, 0.05367636867327377], [-3.020318448930587, -0.657332495932215, 1.2192077454993016], [-3.1314232385001866, -1.5575296086444772, 2.4986162193657244], [-3.378446400292146, -0.8053658014461064, 3.5567896965215997], [-4.105859700541783, -2.4509305976798705, 2.352985539948738], [-1.5468226482012581, -2.509766594543162, 2.821027483166894], [-0.44897476829631316, -1.6306187336043643, 2.6106466879289294], [-1.7496386975713298, -3.255878709799582, 4.0074250275142855], [-1.6320289120976281, -3.5106116624362236, 1.6079357007876873], [-2.345438986660116, 0.44397410980513335, 1.5364663530591296], [-4.247573342377324, -0.3344495623799305, 0.8318812175920258], [-2.596140664371262, -0.7866250332999468, -1.0913507571730208], [-2.6897770261666474, -2.6521285947103417, 0.007528686185012752], [-0.37104516182828323, -1.9490799840121065, 1.2553873021032915], [-0.24552532002048794, -2.1655580448419816, -0.8973245689258732], [1.5770424436171644, 0.0, 0.0], [2.292718146893914, 1.3915527243580579, 0.0], [1.600521547008251, 2.4407219045639, -0.93570860023403], [1.3760692761371693, 1.9145350254105775, -2.1257143081021437], [0.4611709133627504, 2.851946778983252, -0.4188289776134555], [2.4048649076934967, 3.4780800111829113, -1.0689303403306858], [3.540059223330463, 1.2295174124846107, -0.43152105569274557], [2.308846803952248, 1.8960947387583775, 1.2310220414904705], [1.9974224573334811, -0.690678068305522, 1.0535722235493046], [1.9277183224308942, -0.652993231720626, -1.1102241252095317], [-0.3501493572534756, 0.5705349971623068, -1.1530217920585846], [-0.42668432219276003, 0.7576153073313064, 1.0049834283127215], [-1.1175230423143103, -3.176521799292093, 0.8598514408283896]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0474', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
