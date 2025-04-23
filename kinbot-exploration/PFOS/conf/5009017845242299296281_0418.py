import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0418'
logfile = 'conf/5009017845242299296281_0418.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863834, 0.7718203945763811, 1.1635336229088482], [-2.2709622836291934, 0.7431123812655605, 1.1797556627388999], [-2.9970239643019023, 1.465687704195477, 2.3681573973696803], [-2.2704853996782606, 2.7960258391239847, 2.770559335575201], [-3.0948301649158383, 3.5679399740158075, 3.4569311149022734], [-1.1948366276544569, 2.531746623187418, 3.5067904513402466], [-1.692784735145687, 3.7576527233167063, 1.2661660994459563], [-1.4188572202552054, 5.093336324939094, 1.6703245346579063], [-0.7700728935849264, 2.931525692556909, 0.5795266320836723], [-3.0423357848931114, 3.742343803449636, 0.454365489131374], [-4.239914773111565, 1.7601618341371337, 1.9973488528651497], [-3.023054582933852, 0.6702887991292497, 3.4298036615283634], [-2.6427825570536165, 1.3431159399559998, 0.04548369333676551], [-2.6677893778920323, -0.523270637705677, 1.14656640528709], [-0.37104516182828423, 2.06173728721591, 1.0602591291106116], [-0.24552532002048896, 0.3056731502912611, 2.324090564665853], [1.577042443617166, 0.0, 0.0], [2.2927181468939164, 1.3915527243580579, 0.0], [3.782335574419714, 1.3186147352454656, -0.48079377647244215], [4.419592497958776, 0.3492507324843974, 0.14976318527135246], [3.8473329977166397, 1.1104432852892017, -1.7795740852228756], [4.368010411051576, 2.4684340853352094, -0.20560554420872101], [2.293896732220215, 1.8704189044736073, 1.2405689893126324], [1.6494649440008775, 2.235214894314339, -0.8029305726284889], [1.9974224573334833, -0.6906780683055253, 1.0535722235492986], [1.9277183224308958, -0.6529932317206254, -1.1102241252095326], [-0.3501493572534733, -1.283813661620945, 0.08241309473864852], [-0.42668432219275587, 0.4915335256355428, -1.158605816601226], [-3.5877865142349803, 2.989946258831301, 0.7233335604869209]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0418', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
