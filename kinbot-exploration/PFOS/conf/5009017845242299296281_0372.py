import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0372'
logfile = 'conf/5009017845242299296281_0372.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863809, 0.6217394783082152, -1.2501828803164994], [-0.3466020415139029, 2.1278181305643233, -1.566386336981123], [1.1624535490467065, 2.4914860914100156, -1.7945884035232613], [1.9172819807099761, 1.3689994871066287, -2.5882149031527932], [1.146447056654723, 0.9090692451221138, -3.558152857581732], [3.0445827205004443, 1.8471459380411672, -3.1073189910576895], [2.381891830767318, -0.07448949729141746, -1.4827102354710182], [3.5162028704008725, 0.30247148209780267, -0.7122104670133259], [1.1773846443337452, -0.5732038342621308, -0.9297154619172766], [2.8441289587512673, -1.0798258576734945, -2.6035038247709754], [1.2315520144708538, 3.620752460676508, -2.4938554639699135], [1.7685418932081425, 2.6511749791394856, -0.6251499999737744], [-1.011797100425896, 2.4109962356139887, -2.6898422467352834], [-0.8215282967123799, 2.8830315798638466, -0.5833314461381579], [-2.011954187959716, 0.6042709716797812, -1.0485616399675919], [-0.40994706586377744, -0.11587296658230317, -2.3181096973944295], [1.5770424436171684, 0.0, 0.0], [2.2927181468939155, 1.3915527243580579, 0.0], [2.3410798567223194, 2.0598526928949292, 1.4165023767064764], [1.14539164100703, 2.033172555891589, 1.9759511228307933], [3.202217992284646, 1.4422259961598305, 2.1984030628363285], [2.7215555591494938, 3.315205977928797, 1.2745358845394104], [1.629232339193972, 2.21225586731058, -0.8090479336198833], [3.545586830094381, 1.2600392214310825, -0.4280914688619808], [1.9974224573334851, -0.6906780683055236, 1.0535722235493], [1.927718322430895, -0.6529932317206245, -1.1102241252095315], [-0.35014935725347496, 0.7132786644586342, 1.0706086973199325], [-0.4266843221927538, -1.2491488329668505, 0.15362238828850233], [2.602660526153736, -1.9852970993650265, -2.363134914084743]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0372', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
