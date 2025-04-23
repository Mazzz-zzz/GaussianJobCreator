import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1502984803620600000001/kinbot.db')
label = 'conf/1502984803620600000001_r12_insertion_R_6_5_7_0001'
logfile = 'conf/1502984803620600000001_r12_insertion_R_6_5_7_0001.log'

mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(0.0), np.float64(0.0), np.float64(0.0)], [np.float64(-0.35088163747906237), np.float64(0.5207852942357362), np.float64(-1.16847415711209)], [np.float64(-0.5864659549131684), np.float64(0.6913320163695338), np.float64(0.9553547446468357)], [np.float64(-0.42038859015377744), np.float64(-1.2537295953208796), np.float64(0.036496160881157126)], [np.float64(2.015624496583131), np.float64(0.0), np.float64(0.0)], [np.float64(1.8582600963998184), np.float64(-0.48779055842993396), np.float64(1.4377183943705947)], [np.float64(3.4555494140238263), np.float64(-0.3050436261936609), np.float64(0.6161623810060373)], [np.float64(2.0762459987231057), np.float64(1.5981118507965246), np.float64(0.0)], [np.float64(1.395498882618007), np.float64(1.9460329994342869), np.float64(0.5899866267358712)]])

kwargs = {'method': 'bmk', 'basis': '6-31++G(2df,p)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/1502984803620600000001_r12_insertion_R_6_5_7_0001', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999'}
Gaussian.command = 'g16 < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy()  # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    freq = reader_gauss.read_freq(logfile, [np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')])
    zpe = reader_gauss.read_zpe(logfile)
    db.write(mol, name=label, data={'energy': e, 'frequencies': np.asarray(freq),
                                     'zpe': zpe, 'status': 'normal'})

except RuntimeError:
    for i in range(3):
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
            if i == 2:
                db.write(mol, name=label, data={'status': 'error'})
            pass
        else:
            break

with open(logfile, 'a') as f:
    f.write('done\n')
