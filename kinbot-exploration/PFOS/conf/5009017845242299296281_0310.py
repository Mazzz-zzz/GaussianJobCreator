import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0310'
logfile = 'conf/5009017845242299296281_0310.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863828, -1.3935598728845986, 0.08664925740765334], [-0.3466020415139025, -2.420439425248669, -1.059551387211253], [1.1624535490467074, -2.7999021924931133, -1.2603960465750352], [1.924360652536601, -1.7238039150214761, -2.1095693044418966], [1.6581688080331807, -1.8936103403721987, -3.392850108828245], [1.5600180827406158, -0.4992730224820868, -1.7398685800658817], [3.7832317877839654, -1.8440144326786725, -1.8819493028959802], [4.146265721376076, -3.2188116525391695, -1.9058234549119728], [4.364034167085277, -0.8555804182445006, -2.7132275757718536], [3.8645468409417227, -1.3259104534731778, -0.3968571311199683], [1.7462981079693105, -2.8864928243592174, -0.06849875821689058], [1.2535978601361875, -3.964626716020743, -1.8893151990408796], [-1.011797100425896, -3.534969835652363, -0.7430628651027212], [-0.8215282967123789, -1.9466956411138836, -2.2051128650057934], [-2.011954187959718, -1.210216503485697, 0.0009668077395966848], [-0.4099470658637793, -1.9496054034114785, 1.2594037813693577], [1.5770424436171657, 0.0, 0.0], [2.2927181468939177, 1.3915527243580545, 0.0], [2.3410798567223217, 2.0598526928949275, 1.4165023767064737], [1.1453916410070366, 2.0331725558915936, 1.975951122830793], [3.2022179922846448, 1.442225996159824, 2.198403062836329], [2.7215555591495013, 3.3152059779287932, 1.2745358845394037], [1.6292323391939787, 2.2122558673105757, -0.8090479336198896], [3.5455868300943854, 1.2600392214310692, -0.4280914688619795], [1.9974224573334802, -0.6906780683055244, 1.0535722235493032], [1.9277183224308954, -0.6529932317206287, -1.110224125209528], [-0.35014935725347485, 0.5705349971623066, -1.153021792058583], [-0.42668432219275937, 0.7576153073313062, 1.0049834283127228], [4.66816339300337, -0.8031721272569311, -0.2668284179325169]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0310', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
