import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0304'
logfile = 'conf/5009017845242299296281_0304.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863819, 0.7718203945763832, 1.163533622908849], [-2.270962283629193, 0.7431123812655627, 1.1797556627389014], [-2.96999178851271, -0.6604598778102149, 1.2338112746327903], [-2.1956407956954127, -1.6601377768739964, 2.1615018464552516], [-1.781466718272868, -1.034282187843841, 3.2492824118028323], [-2.9823571229144368, -2.675327993379045, 2.506996108856896], [-0.6934923395346745, -2.387095074607015, 1.3031624997854045], [0.1587682255477207, -2.9455546410759648, 2.295233994000793], [-1.1688016120323788, -3.097717967121517, 0.17424577488491666], [-0.054361718048219806, -1.0435727906355243, 0.7864396040672424], [-4.201534988011381, -0.5092377781085995, 1.7127943399246788], [-3.0244357979279344, -1.1831438545952393, 0.015506258743916439], [-2.6212054717929285, 1.4189848017416518, 2.2776416618875825], [-2.7133779807411287, 1.3971249246106954, 0.11259346120335603], [-0.3710451618282825, 2.061737287215913, 1.060259129110612], [-0.24552532002048857, 0.3056731502912613, 2.324090564665855], [1.577042443617165, 0.0, 0.0], [2.2927181468939173, 1.391552724358054, 0.0], [1.6005215470082568, 2.4407219045639, -0.9357086002340298], [1.37606927613717, 1.9145350254105806, -2.1257143081021415], [0.46117091336276284, 2.851946778983257, -0.41882897761345217], [2.4048649076935096, 3.4780800111829127, -1.0689303403306794], [3.5400592233304673, 1.229517412484603, -0.43152105569274385], [2.308846803952255, 1.896094738758371, 1.2310220414904736], [1.997422457333485, -0.6906780683055271, 1.0535722235492977], [1.9277183224308962, -0.6529932317206257, -1.1102241252095313], [-0.3501493572534751, -1.2838136616209437, 0.08241309473865072], [-0.4266843221927575, 0.4915335256355472, -1.1586058166012272], [0.35874395157118266, -1.178567825781499, -0.07789970020524195]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0304', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
