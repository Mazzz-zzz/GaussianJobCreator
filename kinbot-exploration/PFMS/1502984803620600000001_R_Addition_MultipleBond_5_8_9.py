import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFMS/kinbot.db')
label = '1502984803620600000001_R_Addition_MultipleBond_5_8_9'
logfile = '1502984803620600000001_R_Addition_MultipleBond_5_8_9.log'

atom = ['C', 'F', 'F', 'F', 'S', 'O', 'O', 'O', 'H']
geom = [[-0.13873551142583299, -0.10288992553991576, 0.09340764290821346], [-0.6447793155978899, 0.49953254514226586, -0.9666574109735003], [-0.5909044445143531, 0.44023020705502947, 1.1997481009999598], [-0.4080071166535717, -1.387946085362309, 0.06443629241225356], [1.7164679028804517, 0.09626933132876063, 0.04992266142167055], [2.2387224841689504, -0.5549244404318323, 1.1927932380472424], [2.1761926262651863, -0.19949391649403375, -1.2634821660786277], [1.7785881558395047, 1.6640688296420372, 0.17950157866845112], [2.4660002451171934, 1.736047015720779, -1.0209219176982423]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'mp2', 'basis': '6-31G', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_R_Addition_MultipleBond_5_8_9', 'label': '1502984803620600000001_R_Addition_MultipleBond_5_8_9', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'guess': 'Mix,Always', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n5 8 9 F\n'}
Gaussian.command = 'g16 < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy() # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    freq = reader_gauss.read_freq(logfile, ['C', 'F', 'F', 'F', 'S', 'O', 'O', 'O', 'H'])
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
        freq = reader_gauss.read_freq(logfile, ['C', 'F', 'F', 'F', 'S', 'O', 'O', 'O', 'H'])
        zpe = reader_gauss.read_zpe(logfile)
        db.write(mol, name=label, data={'energy': e,
                                         'frequencies': np.asarray(freq),
                                         'zpe': zpe, 'status': 'normal'})
    except RuntimeError:
        db.write(mol, name=label, data={'status': 'error'})

with open(logfile,'a') as f:
    f.write('done\n')
