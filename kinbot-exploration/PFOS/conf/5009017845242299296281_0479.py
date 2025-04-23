import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0479'
logfile = 'conf/5009017845242299296281_0479.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863796, -1.393559872884603, 0.08664925740764864], [-2.2709622836291894, -1.393254564823224, 0.05367636867327149], [-2.9699917885127096, -0.7382819684025511, -1.188880669680407], [-4.434530032795166, -1.266641279561202, -1.3774259648962146], [-4.4121494903032055, -2.456067387613762, -1.9529835970041605], [-5.049857983877611, -1.355022480803869, -0.20164523855452085], [-5.456246996148598, -0.11973325061456502, -2.455553856682981], [-5.840354017548495, 1.0028492669069626, -1.6715249699215797], [-4.791164773781588, -0.026439072178644416, -3.7022770037406376], [-6.711000276427691, -1.0509356712337026, -2.653035215982822], [-3.0171454311974557, 0.5789431357488242, -1.0108747386529344], [-2.2839557356380054, -1.017061933834244, -2.2898193902816026], [-2.621205471792919, -2.6819879407832867, 0.09005605494849198], [-2.7133779807411273, -0.7960712600074827, 1.1536489463716], [-0.3710451618282788, -1.9490799840121193, 1.2553873021032804], [-0.24552532002048516, -2.165558044841981, -0.897324568925883], [1.5770424436171648, 0.0, 0.0], [2.292718146893915, 1.3915527243580559, 0.0], [2.341079856722327, 2.0598526928949403, 1.4165023767064635], [1.1453916410070388, 2.033172555891597, 1.9759511228307907], [3.2022179922846545, 1.4422259961598363, 2.198403062836315], [2.7215555591494924, 3.3152059779287946, 1.2745358845394], [1.6292323391939698, 2.212255867310575, -0.8090479336198902], [3.5455868300943747, 1.26003922143108, -0.428091468861995], [1.9974224573334873, -0.6906780683055228, 1.0535722235492946], [1.9277183224308945, -0.6529932317206173, -1.1102241252095342], [-0.35014935725347596, 0.570534997162309, -1.1530217920585752], [-0.42668432219275915, 0.7576153073312959, 1.0049834283127295], [-6.462610124348825, -1.9786904783885007, -2.536692030331905]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0479', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
