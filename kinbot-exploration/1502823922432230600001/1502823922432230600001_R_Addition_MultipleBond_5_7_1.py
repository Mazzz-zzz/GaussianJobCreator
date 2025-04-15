import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1502823922432230600001/kinbot.db')
label = '1502823922432230600001_R_Addition_MultipleBond_5_7_1'
logfile = '1502823922432230600001_R_Addition_MultipleBond_5_7_1.log'

atom = [np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')]
geom = [[np.float64(-0.014486065476529903), np.float64(-0.11996109275302437), np.float64(1.4437586639259288)], [np.float64(-1.2733173230670713), np.float64(-0.17006732779378536), np.float64(1.77619788399399)], [np.float64(0.5834902225173187), np.float64(0.7718772506261343), np.float64(2.208343005555392)], [np.float64(0.5144507805104077), np.float64(-1.3143159709820258), np.float64(1.6273165832287442)], [np.float64(1.6777428000180659), np.float64(0.17192595558634644), np.float64(0.04771064192393514)], [np.float64(2.0279044074006927), np.float64(-0.6293447030073154), np.float64(-1.085054922711737)], [np.float64(-0.022088841688604872), np.float64(0.255829899279718), np.float64(0.07999750094709183)], [np.float64(1.9299111698059703), np.float64(1.7047493031022087), np.float64(-0.2934688354233593)], [np.float64(1.7723636068630828), np.float64(2.046889883932061), np.float64(-1.1866001016675756)]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'mp2', 'basis': '6-31G', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502823922432230600001_R_Addition_MultipleBond_5_7_1', 'label': '1502823922432230600001_R_Addition_MultipleBond_5_7_1', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'guess': 'Mix,Always', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 7 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n5 7 1 F\n'}
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
