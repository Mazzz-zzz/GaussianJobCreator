import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0288'
logfile = 'conf/5009017845242299296281_0288.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863828, 0.6217394783082179, -1.2501828803165003], [-0.3466020415139025, 2.1278181305643273, -1.56638633698112], [1.1624535490467074, 2.491486091410018, -1.794588403523258], [1.924360652536601, 2.6888425662012905, -0.43807332933071264], [1.6581688080331807, 3.8850995556641563, 0.056510394782913545], [1.5600180827406156, 1.7564069008244552, 0.4375511691392204], [3.7832317877839654, 2.5518231212816724, -0.6559886921968815], [4.146265721376075, 3.259897353351585, -1.8346609336402946], [4.364034167085277, 2.7775142159891413, 0.6156594107056839], [3.8645468409417227, 1.006643583959496, -0.9498435702911299], [1.7462981079693105, 1.5025680769231322, -2.4655267346281255], [1.2535978601361875, 3.618508316135836, -2.488809853075984], [-1.011797100425896, 2.410996235613996, -2.689842246735281], [-0.8215282967123788, 2.883031579863846, -0.5833314461381526], [-2.0119541879597174, 0.6042709716797847, -1.0485616399675892], [-0.4099470658637793, -0.11587296658229951, -2.3181096973944286], [1.5770424436171657, 0.0, 0.0], [2.2927181468939164, 1.3915527243580557, 0.0], [3.7823355744197134, 1.3186147352454594, -0.4807937764724465], [4.4195924979587735, 0.3492507324843891, 0.14976318527134813], [3.847332997716635, 1.1104432852892006, -1.7795740852228823], [4.368010411051579, 2.468434085335203, -0.20560554420872051], [2.293896732220218, 1.870418904473602, 1.2405689893126337], [1.6494649440008815, 2.2352148943143386, -0.8029305726284854], [1.9974224573334842, -0.6906780683055265, 1.053572223549299], [1.9277183224308927, -0.6529932317206248, -1.1102241252095353], [-0.35014935725347357, 0.7132786644586351, 1.0706086973199336], [-0.42668432219275937, -1.2491488329668488, 0.15362238828849986], [4.66816339300337, 0.6326662520096383, -0.5621532568498298]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0288', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
