import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1502984803620600000001/kinbot.db')
label = '1502984803620600000001_r12_insertion_R_8_5_1'
logfile = '1502984803620600000001_r12_insertion_R_8_5_1.log'

atom = [np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')]
geom = [[np.float64(-0.889592), np.float64(-0.117304), np.float64(-0.233938)], [np.float64(-0.943255), np.float64(0.272363), np.float64(-1.522052)], [np.float64(-1.528112), np.float64(0.843783), np.float64(0.476845)], [np.float64(-1.765806), np.float64(-1.149019), np.float64(-0.19329)], [np.float64(0.744951), np.float64(-0.561172), np.float64(0.375374)], [np.float64(0.861057), np.float64(-0.672734), np.float64(1.804527)], [np.float64(1.544833), np.float64(-1.334189), np.float64(-0.51629)], [np.float64(0.620497), np.float64(1.02982), np.float64(-0.116721)], [np.float64(1.345426), np.float64(1.688453), np.float64(-0.083693)]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'b3lyp', 'basis': '6-31+g*', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r12_insertion_R_8_5_1', 'label': '1502984803620600000001_r12_insertion_R_8_5_1', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n'}
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
