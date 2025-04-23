import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0025'
logfile = 'conf/5009017845242299296281_0025.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863845, 0.7718203945763887, 1.163533622908842], [-0.39761971585595823, 2.318865504557597, 1.2514273698287348], [-1.123388212146639, 3.132984570832187, 2.3792791242903557], [-2.585410647330541, 3.522116108662272, 1.96617840498945], [-3.290357209347752, 3.8253563239325916, 3.04209552987404], [-2.5711488993453706, 4.562809768561837, 1.13813652842627], [-3.456801175356166, 2.116729468845753, 1.0789184253851347], [-3.1643747017820463, 0.907046602507038, 1.7670936416270087], [-4.767434488313374, 2.5656004854849113, 0.7854701317520942], [-2.630960974039623, 2.1473609729493885, -0.2618247650338215], [-1.180898235093628, 2.3869269963298314, 3.4788618839288867], [-0.45262953286663876, 4.249846884503501, 2.6305382892183498], [0.9186702629217386, 2.4110951602602593, 1.4604772899922651], [-0.6863393184328448, 2.8685142811653983, 0.07804489435392555], [-0.25604457595343294, 0.2684346684142536, 2.3183978473060165], [-2.0076024771874477, 0.591702034196639, 1.0837604470856887], [1.577042443617166, 0.0, 0.0], [2.2927181468939177, 1.3915527243580557, 0.0], [2.3410798567223234, 2.05985269289493, 1.416502376706472], [1.1453916410070413, 2.0331725558915967, 1.9759511228307893], [3.202217992284643, 1.4422259961598223, 2.19840306283633], [2.7215555591495084, 3.315205977928793, 1.2745358845394015], [1.6292323391939854, 2.212255867310577, -0.8090479336198889], [3.5455868300943854, 1.2600392214310676, -0.42809146886197724], [1.9974224573334802, -0.6906780683055234, 1.0535722235493095], [1.9277183224308971, -0.6529932317206335, -1.1102241252095253], [-0.3501493572534756, -1.2838136616209472, 0.08241309473865782], [-0.42668432219275443, 0.4915335256355364, -1.158605816601231], [-1.9131310404359811, 1.4995794837693166, -0.22924079053342386]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0025', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
Gaussian.command = 'g16 < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy()  # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    freq = reader_gauss.read_freq(logfile, ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'])
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
            freq = reader_gauss.read_freq(logfile, ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'])
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
