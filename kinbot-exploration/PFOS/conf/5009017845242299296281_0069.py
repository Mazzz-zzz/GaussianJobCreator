import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0069'
logfile = 'conf/5009017845242299296281_0069.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863862, -1.393559872884596, 0.08664925740765565], [-2.2709622836291947, -1.393254564823213, 0.05367636867327356], [-2.9970239643019054, -2.783728318379919, 0.08524408716293887], [-3.05161235149174, -3.4433745877869133, -1.3365407299748207], [-4.021172178383425, -2.8960736485529517, -2.048522026827759], [-1.8922874678070294, -3.280072004057724, -1.9675701435284367], [-3.3760113189425525, -5.289049349093577, -1.2369750907842492], [-2.1660027576274805, -5.933555567062759, -0.8589247973664512], [-4.608145274647481, -5.453958580693031, -0.5587364659771561], [-3.634449254937269, -5.557478945820578, -2.767235857663141], [-2.332949719823826, -3.597057518120821, 0.9016340449500815], [-4.241488402883609, -2.631664135309415, 0.5195418500188262], [-2.6427825570536223, -0.7109480038655721, 1.140430677561324], [-2.6677893778920327, -0.7313203152515876, -1.0264488679511479], [-0.37104516182828834, -1.9490799840121042, 1.2553873021032929], [-0.2455253200204909, -2.165558044841981, -0.8973245689258684], [1.5770424436171646, 0.0, 0.0], [2.2927181468939164, 1.3915527243580543, 0.0], [3.782335574419717, 1.3186147352454554, -0.4807937764724423], [4.419592497958776, 0.3492507324843923, 0.14976318527136862], [3.847332997716638, 1.1104432852891937, -1.7795740852228763], [4.368010411051584, 2.4684340853352005, -0.20560554420871335], [2.293896732220209, 1.8704189044736137, 1.2405689893126324], [1.6494649440008793, 2.2352148943143373, -0.8029305726284963], [1.997422457333477, -0.6906780683055216, 1.0535722235493077], [1.9277183224308951, -0.6529932317206313, -1.1102241252095297], [-0.35014935725347124, 0.570534997162306, -1.1530217920585846], [-0.4266843221927622, 0.7576153073313113, 1.0049834283127204], [-2.8228225226921975, -5.861985768069874, -3.1967705237901214]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0069', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
