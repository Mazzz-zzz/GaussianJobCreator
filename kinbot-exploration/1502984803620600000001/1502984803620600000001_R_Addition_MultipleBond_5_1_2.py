import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1502984803620600000001/kinbot.db')
label = '1502984803620600000001_R_Addition_MultipleBond_5_1_2'
logfile = '1502984803620600000001_R_Addition_MultipleBond_5_1_2.log'

atom = [np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')]
geom = [[np.float64(-1.0769379194170108), np.float64(-0.21719535241070784), np.float64(-0.20858047446083813)], [np.float64(-0.6953003743219308), np.float64(-0.035778899002094675), np.float64(-1.557293580504893)], [np.float64(-1.7277852202231596), np.float64(0.8628006399852408), np.float64(0.2159553793643971)], [np.float64(-1.7829622026792613), np.float64(-1.3323858743741455), np.float64(-0.11447108951980389)], [np.float64(0.7304639754246612), np.float64(-0.3585625757920656), np.float64(0.28356019459276305)], [np.float64(0.905614302245621), np.float64(-0.46232716066210994), np.float64(1.7224483215518944)], [np.float64(1.343742935840447), np.float64(-1.3182371084419626), np.float64(-0.6073110960661809)], [np.float64(1.0845116557400403), np.float64(1.1719261254757105), np.float64(-0.15724652637633296)], [np.float64(1.6718620837126725), np.float64(1.6488601461947299), np.float64(0.4623708471265393)]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'mp2', 'basis': '6-31G', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_R_Addition_MultipleBond_5_1_2', 'label': '1502984803620600000001_R_Addition_MultipleBond_5_1_2', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'guess': 'Mix,Always', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n5 1 2 F\n'}
Gaussian.command = 'g16 < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy() # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    freq = reader_gauss.read_freq(logfile, [np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')])
    zpe = reader_gauss.read_zpe(logfile)
    db.write(mol, name=label, data={'energy': e,'frequencies': np.asarray(freq), 'zpe':zpe, 'status': 'normal'})
except RuntimeError:
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
        db.write(mol, name=label, data={'status': 'error'})

with open(logfile,'a') as f:
    f.write('done\n')
