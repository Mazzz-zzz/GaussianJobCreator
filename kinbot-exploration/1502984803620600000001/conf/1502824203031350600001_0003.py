import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1502984803620600000001/kinbot.db')
label = 'conf/1502824203031350600001_0003'
logfile = 'conf/1502824203031350600001_0003.log'

mol = Atoms(symbols=['C', 'F', 'F', 'F', 'S', 'O', 'O', 'O', 'H'], positions=[[np.float64(-0.7710187513442883), np.float64(0.9854618415437302), np.float64(1.5603260594903077)], [np.float64(0.09233233764246429), np.float64(0.9104154444025283), np.float64(2.5442876337362286)], [np.float64(-0.9733868630712574), np.float64(-0.21556533116712728), np.float64(1.0485462901361187)], [np.float64(-1.91561648961404), np.float64(1.4553695221106095), np.float64(2.019986394765196)], [np.float64(-0.19922933980299606), np.float64(2.1809509532371805), np.float64(0.19664790295715096)], [np.float64(0.9663875567193523), np.float64(0.0), np.float64(0.0)], [np.float64(-0.009723003381871237), np.float64(3.422363832948778), np.float64(0.9105828255659292)], [np.float64(1.2490131048205662), np.float64(1.396603990681499), np.float64(0.0)], [np.float64(0.0), np.float64(0.0), np.float64(0.0)]])

kwargs = {'method': 'bmk', 'basis': '6-31++G(2df,p)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/1502824203031350600001_0003', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
