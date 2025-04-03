import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '660901080000000000001_r12_insertion_R_1_2_3'
logfile = '660901080000000000001_r12_insertion_R_1_2_3.log'

atom = [np.str_('F'), np.str_('C'), np.str_('O'), np.str_('F')]
geom = [[np.float64(1.10131), np.float64(-0.434894), np.float64(0.089498)], [np.float64(-0.453287), np.float64(-0.069988), np.float64(-0.399348)], [np.float64(0.453175), np.float64(1.344125), np.float64(0.247656)], [np.float64(-1.111199), np.float64(-0.848479), np.float64(0.450102)]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'b3lyp', 'basis': '6-31+g*', 'nprocshared': 8, 'mem': '700MW', 'chk': '660901080000000000001_r12_insertion_R_1_2_3', 'label': '660901080000000000001_r12_insertion_R_1_2_3', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n2 3 F\n2 4 F\n'}
Gaussian.command = 'g16 < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy() # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    freq = reader_gauss.read_freq(logfile, [np.str_('F'), np.str_('C'), np.str_('O'), np.str_('F')])
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
        freq = reader_gauss.read_freq(logfile, [np.str_('F'), np.str_('C'), np.str_('O'), np.str_('F')])
        zpe = reader_gauss.read_zpe(logfile)
        db.write(mol, name=label, data={'energy': e,
                                         'frequencies': np.asarray(freq),
                                         'zpe': zpe, 'status': 'normal'})
    except RuntimeError:
        db.write(mol, name=label, data={'status': 'error'})

with open(logfile,'a') as f:
    f.write('done\n')
