import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1502824203031350600001/kinbot.db')
label = '1502824203031350600001_intra_OH_migration_Exocyclic_F_5_6_-1'
logfile = '1502824203031350600001_intra_OH_migration_Exocyclic_F_5_6_-1.log'

atom = [np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')]
geom = [[np.float64(2.430067), np.float64(0.904543), np.float64(-2.460064)], [np.float64(2.160069), np.float64(-0.110747), np.float64(-3.264671)], [np.float64(3.33686), np.float64(0.481199), np.float64(-1.59725)], [np.float64(3.048097), np.float64(1.815742), np.float64(-3.212322)], [np.float64(0.920826), np.float64(1.59505), np.float64(-1.617358)], [np.float64(0.142537), np.float64(0.238644), np.float64(-0.950555)], [np.float64(0.179478), np.float64(2.356809), np.float64(-2.506426)], [np.float64(1.141794), np.float64(1.743831), np.float64(-0.092825)], [np.float64(-0.661289), np.float64(0.412199), np.float64(-0.437841)]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'bmk', 'basis': '6-31++G(2df,p)', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502824203031350600001_intra_OH_migration_Exocyclic_F_5_6_-1', 'label': '1502824203031350600001_intra_OH_migration_Exocyclic_F_5_6_-1', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': ''}
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
