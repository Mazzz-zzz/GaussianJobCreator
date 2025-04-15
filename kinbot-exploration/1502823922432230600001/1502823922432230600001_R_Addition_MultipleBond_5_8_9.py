import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1502823922432230600001/kinbot.db')
label = '1502823922432230600001_R_Addition_MultipleBond_5_8_9'
logfile = '1502823922432230600001_R_Addition_MultipleBond_5_8_9.log'

atom = [np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')]
geom = [[np.float64(-0.7091651562830197), np.float64(-0.09028667905935793), np.float64(1.2470157922439928)], [np.float64(-1.98233643650738), np.float64(-0.0015670003801782286), np.float64(0.9313027444618629)], [np.float64(-0.4182259122282983), np.float64(0.8430857638918406), np.float64(2.146291492382258)], [np.float64(-0.4869704367874309), np.float64(-1.269358757309618), np.float64(1.8176271196323488)], [np.float64(1.7306535414882258), np.float64(0.05153688681792921), np.float64(0.06792691732881605)], [np.float64(1.9944064391757124), np.float64(-0.391467977253563), np.float64(-1.2702598374570855)], [np.float64(0.03002247241701304), np.float64(0.05690092055860781), np.float64(0.12551274166786486)], [np.float64(1.9254100617452987), np.float64(1.6277064296644663), np.float64(0.0012586115882306252)], [np.float64(1.7637775179632853), np.float64(1.5895002929770266), np.float64(-1.3741826834203323)]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'mp2', 'basis': '6-31G', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502823922432230600001_R_Addition_MultipleBond_5_8_9', 'label': '1502823922432230600001_R_Addition_MultipleBond_5_8_9', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'guess': 'Mix,Always', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 7 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n5 8 9 F\n'}
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
