import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '660901080000000000001_well_mp2'
logfile = '660901080000000000001_well_mp2.log'

mol = Atoms(symbols=[np.str_('F'), np.str_('C'), np.str_('O'), np.str_('F')], positions=[[np.float64(1.066168), np.float64(-0.691001), np.float64(0.139573)], [np.float64(-0.002178), np.float64(0.094228), np.float64(0.09028)], [np.float64(0.001777), np.float64(1.270752), np.float64(0.008996)], [np.float64(-1.075767), np.float64(-0.683214), np.float64(0.149058)]])

kwargs = {'method': 'mp2', 'basis': '6-31G', 'nprocshared': 8, 'mem': '700MW', 'chk': '660901080000000000001_well_mp2', 'label': '660901080000000000001_well_mp2', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
Gaussian.command = 'g16 < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy()  # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    freq = reader_gauss.read_freq(logfile, [np.str_('F'), np.str_('C'), np.str_('O'), np.str_('F')])
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
            freq = reader_gauss.read_freq(logfile, [np.str_('F'), np.str_('C'), np.str_('O'), np.str_('F')])
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
