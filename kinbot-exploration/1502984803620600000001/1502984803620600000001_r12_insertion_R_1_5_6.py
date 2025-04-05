import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1502984803620600000001/kinbot.db')
label = '1502984803620600000001_r12_insertion_R_1_5_6'
logfile = '1502984803620600000001_r12_insertion_R_1_5_6.log'

atom = [np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')]
geom = [[np.float64(-0.794629), np.float64(-0.22338), np.float64(-0.155724)], [np.float64(-1.261648), np.float64(-0.040658), np.float64(-1.449248)], [np.float64(-1.28305), np.float64(0.92382), np.float64(0.408272)], [np.float64(-1.605083), np.float64(-1.239287), np.float64(0.229869)], [np.float64(0.855602), np.float64(-0.430778), np.float64(-0.005306)], [np.float64(-0.026169), np.float64(-0.477741), np.float64(1.563218)], [np.float64(1.559215), np.float64(-1.271532), np.float64(-0.827549)], [np.float64(1.469543), np.float64(1.018749), np.float64(-0.151892)], [np.float64(1.076222), np.float64(1.740804), np.float64(0.379126)]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'b3lyp', 'basis': '6-31+g*', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r12_insertion_R_1_5_6', 'label': '1502984803620600000001_r12_insertion_R_1_5_6', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n'}
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
