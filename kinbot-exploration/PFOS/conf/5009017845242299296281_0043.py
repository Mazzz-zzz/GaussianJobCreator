import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0043'
logfile = 'conf/5009017845242299296281_0043.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863792, 0.6217394783082055, -1.250182880316509], [-2.2709622836291903, 0.6501421835576418, -1.2334320314121892], [-2.9699917885127127, 1.3987418462127552, -0.04493060495240272], [-4.434530032795169, 1.8262065172130013, -0.4082305431339445], [-5.030583914311103, 0.8584013505576134, -1.0824522224071096], [-5.127712885957738, 2.088652764913049, 0.6960880430098283], [-4.466501034862393, 3.375387929270144, -1.4668514094871647], [-5.74752992127994, 3.470779758288445, -2.0770123709804125], [-3.873416736570148, 4.408723016758298, -0.701458449150409], [-3.4203658901754115, 2.928287179206995, -2.5560000715864777], [-3.0171454311974597, 0.5859716358429843, 1.006816832231561], [-2.283955735638011, 2.491572728979181, 0.2641082232181906], [-2.6212054717929196, 1.2630031390416094, -2.3676977168361044], [-2.7133779807411247, -0.6010536646032352, -1.2662424075749648], [-0.3710451618282757, -0.1126573032038186, -2.3156464312139025], [-0.24552532002048638, 1.8598848945507076, -1.4267659957399927], [1.5770424436171637, 0.0, 0.0], [2.292718146893914, 1.391552724358052, 0.0], [3.782335574419714, 1.3186147352454543, -0.48079377647244725], [4.419592497958774, 0.34925073248439753, 0.1497631852713588], [3.8473329977166433, 1.1104432852891855, -1.7795740852228739], [4.368010411051579, 2.468434085335205, -0.20560554420872296], [2.2938967322202104, 1.8704189044736128, 1.2405689893126257], [1.6494649440008802, 2.2352148943143337, -0.8029305726284972], [1.9974224573334818, -0.6906780683055209, 1.0535722235493046], [1.9277183224308994, -0.6529932317206317, -1.1102241252095226], [-0.35014935725348, 0.7132786644586374, 1.070608697319925], [-0.42668432219275776, -1.2491488329668508, 0.1536223882885046], [-3.8742795550373583, 2.558007305015701, -3.325946554347763]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0043', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
