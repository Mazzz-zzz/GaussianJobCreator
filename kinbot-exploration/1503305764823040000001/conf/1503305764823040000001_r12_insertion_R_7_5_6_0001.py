import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/1503305764823040000001/kinbot.db')
label = 'conf/1503305764823040000001_r12_insertion_R_7_5_6_0001'
logfile = 'conf/1503305764823040000001_r12_insertion_R_7_5_6_0001.log'

mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(0.0), np.float64(0.0), np.float64(0.0)], [np.float64(-0.34684157427781814), np.float64(0.5096695864407945), np.float64(-1.1747377811915107)], [np.float64(-0.5882457192333358), np.float64(0.7008306830462877), np.float64(0.9469056619640256)], [np.float64(-0.4211046145682315), np.float64(-1.2529448339110907), np.float64(0.0471213594652229)], [np.float64(2.0131464862649713), np.float64(0.0), np.float64(0.0)], [np.float64(1.8670038807160683), np.float64(-0.4875885231064331), np.float64(1.4380525491261666)], [np.float64(3.459662395968114), np.float64(-0.3042700748564122), np.float64(0.6051907235678371)], [np.float64(2.072940182433764), np.float64(1.598077028089533), np.float64(0.0)], [np.float64(1.3958175688913315), np.float64(1.9450786587506723), np.float64(0.5947221793879692)]])

kwargs = {'method': 'bmk', 'basis': '6-31++G(2df,p)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/1503305764823040000001_r12_insertion_R_7_5_6_0001', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999'}
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
