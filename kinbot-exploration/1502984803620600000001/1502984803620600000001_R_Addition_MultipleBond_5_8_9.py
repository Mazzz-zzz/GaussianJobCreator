import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1502984803620600000001/kinbot.db')
label = '1502984803620600000001_R_Addition_MultipleBond_5_8_9'
logfile = '1502984803620600000001_R_Addition_MultipleBond_5_8_9.log'

atom = [np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')]
geom = [[np.float64(-1.1161243257488833), np.float64(-0.21124870576759489), np.float64(-0.33634928621579524)], [np.float64(-1.1657034768813903), np.float64(0.03669043004598813), np.float64(-1.6439444327240058)], [np.float64(-1.679665434947475), np.float64(0.7847870238565526), np.float64(0.36008520955275414)], [np.float64(-1.7124662470001757), np.float64(-1.3706305128203449), np.float64(-0.055766249402511306)], [np.float64(0.6917122830820852), np.float64(-0.34755880528537303), np.float64(0.18174659262722034)], [np.float64(0.7767106261014463), np.float64(-0.2777046780075689), np.float64(1.6330671302673168)], [np.float64(1.2833327231361895), np.float64(-1.4312996141623144), np.float64(-0.5722408296592154)], [np.float64(1.1768977234991307), np.float64(1.1134074429940362), np.float64(-0.3554729887972055)], [np.float64(1.6166102820735286), np.float64(1.435927892214598), np.float64(0.9187427057702715)]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'mp2', 'basis': '6-31G', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_R_Addition_MultipleBond_5_8_9', 'label': '1502984803620600000001_R_Addition_MultipleBond_5_8_9', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'guess': 'Mix,Always', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n5 8 9 F\n'}
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
