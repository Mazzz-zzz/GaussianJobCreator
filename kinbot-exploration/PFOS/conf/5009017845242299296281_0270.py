import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0270'
logfile = 'conf/5009017845242299296281_0270.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863809, -1.3935598728846008, 0.08664925740765222], [-0.3466020415139029, -2.4204394252486647, -1.0595513872112599], [-0.7363023803695192, -2.0273092602547944, -2.527598274606109], [-0.5090926255515552, -0.49965098445894157, -2.800033732485606], [-0.4446952607212935, -0.2804007433113329, -4.101680479585844], [-1.5014279159495576, 0.21624569992661052, -2.2788210606620134], [1.0867274347625948, 0.11649968183315912, -2.028425008687243], [0.8713180636927215, 0.28072126617334536, -0.6323372313959743], [2.1365540635068037, -0.6546895103698991, -2.583842015139847], [1.1313842640058616, 1.5388540618870508, -2.703300119320648], [0.018370741190132885, -2.7227744393129885, -3.373640954864597], [-2.0169579203546344, -2.300385488260373, -2.741562204025709], [0.9809271994195095, -2.56232718372628, -1.010956477982447], [-0.9145558447621676, -3.5839764719224814, -0.7664486088494608], [-2.0119541879597143, -1.2102165034857024, 0.0009668077396009885], [-0.4099470658637716, -1.9496054034114867, 1.2594037813693533], [1.577042443617166, 0.0, 0.0], [2.292718146893914, 1.391552724358057, 0.0], [1.6005215470082597, 2.4407219045638953, -0.9357086002340256], [1.376069276137155, 1.9145350254105828, -2.125714308102139], [0.4611709133627484, 2.851946778983252, -0.41882897761344584], [2.4048649076934883, 3.4780800111829095, -1.0689303403306907], [3.5400592233304646, 1.2295174124846093, -0.4315210556927474], [2.3088468039522545, 1.896094738758375, 1.2310220414904698], [1.9974224573334864, -0.6906780683055288, 1.0535722235492953], [1.9277183224308947, -0.6529932317206225, -1.1102241252095344], [-0.35014935725347784, 0.5705349971623097, -1.1530217920585788], [-0.4266843221927562, 0.7576153073313006, 1.0049834283127304], [2.045033102647419, 1.789061016861257, -2.899780820732956]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0270', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
