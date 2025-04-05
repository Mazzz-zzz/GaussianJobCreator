import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1502984803620600000001/kinbot.db')
label = '1502984803620600000001_well'
logfile = '1502984803620600000001_well.log'

mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(-1.02), np.float64(-0.175481808484), np.float64(-0.286312424369)], [np.float64(-1.1), np.float64(0.009235884657), np.float64(-1.597808045674)], [np.float64(-1.54), np.float64(0.877409042422), np.float64(0.332491847655)], [np.float64(-1.68), np.float64(-1.256080313362), np.float64(0.046179423285)], [np.float64(0.74), np.float64(-0.369435386283), np.float64(0.212425347113)], [np.float64(0.75), np.float64(-0.452558348197), np.float64(1.625515699645)], [np.float64(1.31), np.float64(-1.339203275276), np.float64(-0.637276041338)], [np.float64(1.29), np.float64(1.043654966249), np.float64(-0.22166123177)], [np.float64(1.24), np.float64(1.662459238273), np.float64(0.517209540796)]])

kwargs = {'method': 'b3lyp', 'basis': '6-31+g*', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_well', 'label': '1502984803620600000001_well', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
Gaussian.command = 'g16 < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy()  # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    freq = reader_gauss.read_freq(logfile, [np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')])
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
            freq = reader_gauss.read_freq(logfile, [np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')])
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
