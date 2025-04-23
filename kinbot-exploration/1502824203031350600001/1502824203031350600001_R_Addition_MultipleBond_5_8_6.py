import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1502824203031350600001/kinbot.db')
label = '1502824203031350600001_R_Addition_MultipleBond_5_8_6'
logfile = '1502824203031350600001_R_Addition_MultipleBond_5_8_6.log'

atom = [np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')]
geom = [[np.float64(2.2468108041802024), np.float64(0.8753841263429437), np.float64(-2.466860870582462)], [np.float64(1.6141686208713868), np.float64(-0.10458661326575139), np.float64(-3.0596856134052337)], [np.float64(3.1468896587411903), np.float64(0.4295625129081621), np.float64(-1.6180727261456107)], [np.float64(2.8558817822426033), np.float64(1.61495733774928), np.float64(-3.384026480984064)], [np.float64(1.0538688417260265), np.float64(2.1069827651191426), np.float64(-1.62353948745526)], [np.float64(0.6521245546512807), np.float64(0.24476119965719897), np.float64(-0.598632325316686)], [np.float64(-0.2351664589280126), np.float64(1.928911962183809), np.float64(-2.267436014732452)], [np.float64(1.1283689182438275), np.float64(1.4920178670974982), np.float64(-0.11401645016883936)], [np.float64(0.2731155118131982), np.float64(-0.226923385804516), np.float64(0.15840469472833968)]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'mp2', 'basis': '6-31G', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502824203031350600001_R_Addition_MultipleBond_5_8_6', 'label': '1502824203031350600001_R_Addition_MultipleBond_5_8_6', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'guess': 'Mix,Always', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 7 F\n5 8 F\n6 8 F\n6 9 F\n5 8 6 F\n'}
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
