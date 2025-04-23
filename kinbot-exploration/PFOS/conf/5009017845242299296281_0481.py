import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0481'
logfile = 'conf/5009017845242299296281_0481.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586383, 0.6217394783082151, -1.2501828803165005], [-0.39761971585595646, -0.07566485901595643, -2.633910119820667], [1.087941409756326, -0.0836723585219196, -3.1385123060458486], [1.4730420899455359, 1.2684473752977057, -3.8335337809193666], [0.9614495666061597, 2.28235614395524, -3.1576905727862896], [2.7953494776680783, 1.3981654073395335, -3.8933371678685162], [0.8248630198382829, 1.357801655147954, -5.592380856950236], [-0.5063280479077278, 0.8575285521396571, -5.601770818061772], [1.20042992580093, 2.623237073570472, -6.10519609911399], [1.7738962514906855, 0.28030170799658344, -6.239584394422953], [1.8947536492841413, -0.26661420225360816, -2.0970816341016434], [1.259810370563856, -1.0641590715978366, -4.015744696739174], [-0.7843720924895171, -1.3445493546988774, -2.475221044816154], [-1.1581983763610537, 0.4930222734242938, -3.561495300206357], [-0.2560445759534288, 1.8735740976390431, -1.3916701657561934], [-2.007602477187445, 0.6427130616946763, -1.0543092166280599], [1.5770424436171646, 0.0, 0.0], [2.292718146893914, 1.3915527243580577, 0.0], [2.3410798567223203, 2.0598526928949337, 1.4165023767064788], [1.1453916410070346, 2.033172555891593, 1.9759511228307973], [3.2022179922846425, 1.4422259961598294, 2.1984030628363254], [2.721555559149497, 3.315205977928801, 1.2745358845394015], [1.6292323391939751, 2.2122558673105797, -0.8090479336198894], [3.5455868300943787, 1.260039221431079, -0.4280914688619839], [1.9974224573334838, -0.6906780683055236, 1.0535722235492992], [1.9277183224308958, -0.6529932317206231, -1.1102241252095322], [-0.35014935725347696, 0.7132786644586339, 1.0706086973199318], [-0.42668432219275315, -1.2491488329668503, 0.15362238828850344], [1.3346905127083428, -0.5812711051179595, -6.267025119970061]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0481', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
