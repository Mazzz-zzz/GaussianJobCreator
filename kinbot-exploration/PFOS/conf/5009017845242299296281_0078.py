import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0078'
logfile = 'conf/5009017845242299296281_0078.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863834, 0.6217394783082084, -1.2501828803165063], [-2.2709622836291925, 0.6501421835576441, -1.233432031412183], [-2.969991788512714, 1.3987418462127543, -0.04493060495239275], [-2.1956407956954216, 2.7019843977942113, 0.3569705653274638], [-2.991273105863037, 3.5033560351048902, 1.04346504382554], [-1.1300133907845076, 2.3990794081170197, 1.0928964051269026], [-1.584021457030936, 3.641795456090555, -1.147801939202883], [-0.43246898602355366, 2.9756849415392463, -1.6504496256901326], [-2.723451699091648, 3.9546563740467544, -1.9282976413954638], [-1.126754123236767, 4.962120354252738, -0.42103417812181243], [-4.201534988011384, 1.73794229888726, -0.4153843175535624], [-3.0244357979279393, 0.6050007412875007, 1.0168795050389519], [-2.621205471792928, 1.2630031390416168, -2.3676977168360955], [-2.7133779807411287, -0.6010536646032322, -1.2662424075749599], [-0.3710451618282809, -0.11265730320380932, -2.3156464312139025], [-0.2455253200204912, 1.859884894550714, -1.426765995739986], [1.5770424436171644, 0.0, 0.0], [2.292718146893913, 1.3915527243580588, 0.0], [3.7823355744197116, 1.318614735245462, -0.48079377647244753], [4.419592497958774, 0.3492507324844, 0.1497631852713564], [3.8473329977166366, 1.110443285289201, -1.7795740852228763], [4.368010411051576, 2.4684340853352085, -0.20560554420871996], [2.2938967322202104, 1.8704189044736113, 1.2405689893126286], [1.6494649440008766, 2.2352148943143364, -0.8029305726284909], [1.997422457333483, -0.6906780683055241, 1.0535722235493021], [1.927718322430899, -0.6529932317206246, -1.1102241252095266], [-0.3501493572534803, 0.7132786644586355, 1.0706086973199305], [-0.4266843221927585, -1.2491488329668534, 0.15362238828850155], [-1.608225876241212, 5.071806319019493, 0.4109011119627079]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0078', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
