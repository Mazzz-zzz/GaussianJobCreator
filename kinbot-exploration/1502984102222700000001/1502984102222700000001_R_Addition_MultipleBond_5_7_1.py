import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1502984102222700000001/kinbot.db')
label = '1502984102222700000001_R_Addition_MultipleBond_5_7_1'
logfile = '1502984102222700000001_R_Addition_MultipleBond_5_7_1.log'

atom = [np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')]
geom = [[np.float64(1.130592141493073), np.float64(1.6932674719497005), np.float64(1.2739312611327953)], [np.float64(-0.13107835497749032), np.float64(1.602261585896308), np.float64(0.8947047554098069)], [np.float64(1.3333278044095573), np.float64(2.92551311554895), np.float64(1.639213538354292)], [np.float64(1.3573895404950227), np.float64(0.9357720108217472), np.float64(2.312898477358915)], [np.float64(1.5663298870663962), np.float64(-0.04271190138376623), np.float64(0.08135031604283052)], [np.float64(2.6445253113116434), np.float64(-0.8112955886035931), np.float64(0.5926916302901043)], [np.float64(2.0344301674360437), np.float64(1.502334591442027), np.float64(0.20033819272229017)], [np.float64(1.0713236463940854), np.float64(-0.17595665702185065), np.float64(-1.2381916011447116)], [np.float64(0.5077514858219063), np.float64(-0.28862989639151604), np.float64(0.8905427455627194)]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'mp2', 'basis': '6-31G', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984102222700000001_R_Addition_MultipleBond_5_7_1', 'label': '1502984102222700000001_R_Addition_MultipleBond_5_7_1', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'guess': 'Mix,Always', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 7 F\n5 6 F\n5 7 F\n5 8 F\n5 9 F\n5 7 1 F\n'}
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
