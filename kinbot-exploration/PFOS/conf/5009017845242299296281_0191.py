import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0191'
logfile = 'conf/5009017845242299296281_0191.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863796, -1.393559872884603, 0.08664925740764864], [-2.2709622836291894, -1.393254564823224, 0.05367636867327149], [-2.9699917885127096, -0.7382819684025511, -1.188880669680407], [-4.434530032795166, -1.266641279561202, -1.3774259648962146], [-4.4121494903032055, -2.456067387613762, -1.9529835970041605], [-5.049857983877611, -1.355022480803869, -0.20164523855452085], [-5.456246996148598, -0.11973325061456502, -2.455553856682981], [-5.840354017548494, 1.0028492669069633, -1.6715249699215793], [-4.791164773781588, -0.026439072178644416, -3.7022770037406376], [-6.711000276427691, -1.0509356712337026, -2.653035215982822], [-3.0171454311974557, 0.5789431357488242, -1.0108747386529344], [-2.2839557356380054, -1.017061933834244, -2.2898193902816026], [-2.621205471792919, -2.6819879407832867, 0.09005605494849198], [-2.7133779807411273, -0.7960712600074827, 1.1536489463716], [-0.3710451618282788, -1.9490799840121193, 1.2553873021032804], [-0.24552532002048516, -2.165558044841981, -0.897324568925883], [1.5770424436171648, 0.0, 0.0], [2.292718146893915, 1.3915527243580559, 0.0], [3.7823355744197116, 1.318614735245476, -0.4807937764724385], [4.4195924979587735, 0.349250732484399, 0.1497631852713519], [3.847332997716631, 1.1104432852892145, -1.7795740852228752], [4.368010411051576, 2.468434085335202, -0.20560554420871074], [2.2938967322202135, 1.870418904473604, 1.2405689893126375], [1.6494649440008744, 2.235214894314342, -0.8029305726284759], [1.9974224573334862, -0.6906780683055224, 1.053572223549295], [1.9277183224308947, -0.6529932317206179, -1.110224125209534], [-0.35014935725347596, 0.570534997162309, -1.1530217920585752], [-0.42668432219275915, 0.7576153073312959, 1.0049834283127295], [-7.076477711140584, -0.9383983122700577, -3.5416993675323107]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0191', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
