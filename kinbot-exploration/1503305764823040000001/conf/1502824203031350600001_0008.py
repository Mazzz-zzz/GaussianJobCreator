import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1503305764823040000001/kinbot.db')
label = 'conf/1502824203031350600001_0008'
logfile = 'conf/1502824203031350600001_0008.log'

mol = Atoms(symbols=['C', 'F', 'F', 'F', 'S', 'O', 'O', 'O', 'H'], positions=[[np.float64(-0.09700088620247037), np.float64(-1.86639534467334), np.float64(-0.3555905615446222)], [np.float64(-1.376841177999895), np.float64(-2.1877075992730193), np.float64(-0.3036656500477846)], [np.float64(0.5608109522768345), np.float64(-2.567205973910476), np.float64(0.5372128150460612)], [np.float64(0.36542250090409933), np.float64(-2.126300690643457), np.float64(-1.5619755634055705)], [np.float64(0.0), np.float64(0.0), np.float64(0.0)], [np.float64(1.656630526535413), np.float64(0.0), np.float64(0.0)], [np.float64(2.056917267801028), np.float64(1.3657683715451001), np.float64(0.0)], [np.float64(-0.4877078568772406), np.float64(0.10816016320754973), np.float64(1.3592301325795146)], [np.float64(2.2947747443043), np.float64(1.4995218233139247), np.float64(0.9280048747918408)]])

kwargs = {'method': 'bmk', 'basis': '6-31++G(2df,p)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/1502824203031350600001_0008', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
Gaussian.command = 'g16 < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy()  # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    freq = reader_gauss.read_freq(logfile, ['C', 'F', 'F', 'F', 'S', 'O', 'O', 'O', 'H'])
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
            freq = reader_gauss.read_freq(logfile, ['C', 'F', 'F', 'F', 'S', 'O', 'O', 'O', 'H'])
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
