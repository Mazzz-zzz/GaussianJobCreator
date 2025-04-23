import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0384'
logfile = 'conf/5009017845242299296281_0384.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863841, 0.6217394783082171, -1.2501828803165], [-0.34660204151390467, 2.1278181305643282, -1.5663863369811164], [1.1624535490467052, 2.4914860914100205, -1.7945884035232549], [1.3229555267971602, 3.820814416129085, -2.611111813622318], [0.43579710458927495, 4.710009286593795, -2.2002780767872596], [2.5484901422117674, 4.313914635838239, -2.456783850192621], [1.0607552256216175, 3.5463180474416416, -4.44892738564141], [0.8244098215666444, 4.809473170583544, -5.058034102623292], [2.0568619372694896, 2.6342438976813116, -4.8746607951796435], [-0.30216120032047983, 2.761022150210085, -4.369480551233039], [1.750051457246917, 2.6500962336579237, -0.6119798150372605], [1.7633282149005813, 1.5176064249120498, -2.4661074521032558], [-1.011797100425899, 2.4109962356139967, -2.6898422467352776], [-0.8215282967123796, 2.883031579863847, -0.5833314461381487], [-2.011954187959719, 0.604270971679783, -1.0485616399675886], [-0.40994706586378016, -0.115872966582296, -2.3181096973944304], [1.5770424436171662, 0.0, 0.0], [2.2927181468939164, 1.391552724358055, 0.0], [2.3410798567223274, 2.0598526928949266, 1.416502376706474], [1.1453916410070404, 2.03317255589159, 1.9759511228307949], [3.2022179922846523, 1.4422259961598285, 2.1984030628363227], [2.7215555591495013, 3.315205977928795, 1.274535884539403], [1.6292323391939707, 2.212255867310574, -0.809047933619886], [3.5455868300943796, 1.2600392214310763, -0.42809146886198546], [1.9974224573334833, -0.6906780683055298, 1.0535722235492964], [1.927718322430893, -0.6529932317206233, -1.1102241252095377], [-0.35014935725347385, 0.7132786644586336, 1.0706086973199358], [-0.4266843221927559, -1.249148832966851, 0.15362238828849786], [-0.39305727030787, 2.332027048358138, -3.5071212909801663]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0384', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
