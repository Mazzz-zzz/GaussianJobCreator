import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0280'
logfile = 'conf/5009017845242299296281_0280.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586382, 0.621739478308214, -1.2501828803165023], [-0.3466020415139027, 2.127818130564323, -1.5663863369811213], [1.1624535490467065, 2.491486091410016, -1.79458840352326], [1.9243606525365986, 2.688842566201291, -0.438073329330714], [1.5424767514304087, 1.7679870357247436, 0.429504748506474], [3.2378450948526947, 2.601143892671068, -0.6277245961538497], [1.5819802655580875, 4.373096286345696, 0.31540421373688055], [2.3511039545279884, 5.340046772014328, -0.3887878887042926], [0.1817591264595308, 4.460393050776604, 0.5080614435245339], [2.2429113162508125, 4.130640923929035, 1.7242790035158122], [1.7462981079693156, 1.5025680769231307, -2.465526734628126], [1.2535978601361852, 3.618508316135834, -2.488809853075987], [-1.011797100425897, 2.4109962356139887, -2.6898422467352834], [-0.8215282967123817, 2.883031579863842, -0.5833314461381585], [-2.0119541879597156, 0.6042709716797781, -1.0485616399675917], [-0.40994706586377727, -0.11587296658230416, -2.3181096973944286], [1.577042443617168, 0.0, 0.0], [2.2927181468939155, 1.3915527243580572, 0.0], [2.3410798567223208, 2.059852692894929, 1.4165023767064742], [1.1453916410070328, 2.0331725558915896, 1.9759511228307935], [3.2022179922846465, 1.442225996159829, 2.1984030628363276], [2.7215555591494933, 3.315205977928797, 1.2745358845394066], [1.6292323391939711, 2.2122558673105743, -0.8090479336198866], [3.5455868300943796, 1.260039221431081, -0.4280914688619823], [1.997422457333486, -0.6906780683055258, 1.0535722235492988], [1.9277183224308954, -0.6529932317206254, -1.1102241252095324], [-0.3501493572534749, 0.7132786644586344, 1.0706086973199327], [-0.4266843221927563, -1.2491488329668516, 0.15362238828850022], [3.1510582831678966, 4.464080638833583, 1.731309629353699]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0280', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
