import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1502824203031350600001/kinbot.db')
label = '1502824203031350600001_R_Addition_MultipleBond_5_1_2'
logfile = '1502824203031350600001_R_Addition_MultipleBond_5_1_2.log'

atom = [np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')]
geom = [[np.float64(2.3068057543185625), np.float64(1.044893915489242), np.float64(-2.403744602096733)], [np.float64(1.314238462871646), np.float64(0.10750895584912563), np.float64(-2.7622478368888626)], [np.float64(3.2198212946530598), np.float64(0.42989158344053846), np.float64(-1.7031456454450757)], [np.float64(2.765737952040744), np.float64(1.5613132932056653), np.float64(-3.521195356245281)], [np.float64(0.9978826678707916), np.float64(2.087238985916737), np.float64(-1.506775067428973)], [np.float64(0.9419306463982363), np.float64(-0.038448106763306844), np.float64(-0.0034175431647686435)], [np.float64(-0.27079218222752116), np.float64(1.7802532125991042), np.float64(-2.1428383096168786)], [np.float64(1.1953910519484043), np.float64(1.3657092493204055), np.float64(-0.05512148556789301)], [np.float64(0.08823020651757942), np.float64(-0.09731353435232008), np.float64(0.4529686039613489)]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'mp2', 'basis': '6-31G', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502824203031350600001_R_Addition_MultipleBond_5_1_2', 'label': '1502824203031350600001_R_Addition_MultipleBond_5_1_2', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'guess': 'Mix,Always', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 7 F\n5 8 F\n6 8 F\n6 9 F\n5 1 2 F\n'}
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
