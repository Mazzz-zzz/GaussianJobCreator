import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0154'
logfile = 'conf/5009017845242299296281_0154.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586384, 0.621739478308218, -1.2501828803164992], [-0.39761971585595884, -0.07566485901595232, -2.6339101198206665], [-1.123388212146639, 0.4940238789133734, -3.9028837901505353], [-2.5854106473305416, -0.0582976072378755, -4.033331227674647], [-3.290357209347749, 0.7218538476437423, -4.833903519990093], [-2.5711488993453706, -1.2957497376887053, -4.5205774364234905], [-3.456801175356166, -0.1239939694282254, -2.372600705652132], [-4.854442543346009, -0.24073142942566755, -2.6079769570600067], [-2.713755776739175, -1.0080908700177664, -1.5530807378599571], [-3.14607173426834, 1.3441120826462953, -1.8944224291564145], [-1.1808982350936263, 1.8193192695749194, -3.8065703577649503], [-0.45262953286663765, 0.15318954183903086, -4.995744508783358], [0.9186702629217389, 0.05926285465344136, -2.8183083047232302], [-0.686339318432846, -1.3666682794365106, -2.523228685784664], [-0.2560445759534301, 1.873574097639047, -1.391670165756191], [-2.007602477187447, 0.6427130616946776, -1.0543092166280577], [1.5770424436171657, 0.0, 0.0], [2.2927181468939146, 1.3915527243580592, 0.0], [3.7823355744197116, 1.318614735245466, -0.4807937764724492], [4.419592497958777, 0.3492507324843964, 0.1497631852713472], [3.847332997716634, 1.1104432852892068, -1.7795740852228805], [4.368010411051577, 2.468434085335211, -0.20560554420871857], [2.29389673222022, 1.870418904473607, 1.2405689893126326], [1.6494649440008802, 2.2352148943143426, -0.8029305726284883], [1.9974224573334864, -0.6906780683055279, 1.0535722235492966], [1.927718322430894, -0.6529932317206213, -1.1102241252095353], [-0.3501493572534722, 0.7132786644586361, 1.0706086973199334], [-0.42668432219275626, -1.2491488329668494, 0.1536223882885022], [-3.8863597961933274, 1.9300833995242275, -2.105547008170058]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0154', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
