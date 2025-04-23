import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0028'
logfile = 'conf/5009017845242299296281_0028.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863845, 0.6217394783082163, -1.2501828803165014], [-0.39761971585595723, -0.07566485901595355, -2.6339101198206682], [-0.7491833517666582, -1.5988198371506652, -2.7680043658220987], [0.3751315092500395, -2.510239815805813, -2.1640052972737625], [-0.10815375808507399, -3.715131603382632, -1.91667338731517], [1.3965384052652348, -2.60780306994624, -3.0102377740596302], [1.0485429431935351, -1.8142296954322228, -0.5566004683377848], [1.9279419990532713, -0.7415155439916472, -0.8702446260422215], [-0.04579798938659392, -1.7078876223035397, 0.33584900273066], [1.907047804188475, -3.0576466343715816, -0.11227361857993147], [-1.8820373528306367, -1.8378201814607684, -2.1133959330938117], [-0.8963077255402105, -1.9155699118782812, -4.048025307170635], [-1.1305489784219283, 0.5952337900286083, -3.5269594944504696], [0.8878652548597807, 0.09037708117406702, -2.9207161247297284], [-0.2560445759534276, 1.8735740976390458, -1.3916701657561938], [-2.007602477187446, 0.6427130616946759, -1.0543092166280634], [1.577042443617165, 0.0, 0.0], [2.2927181468939137, 1.3915527243580577, 0.0], [3.782335574419713, 1.3186147352454616, -0.48079377647244353], [4.4195924979587735, 0.34925073248439464, 0.14976318527135052], [3.847332997716633, 1.1104432852892039, -1.7795740852228754], [4.368010411051577, 2.468434085335206, -0.20560554420870952], [2.2938967322202135, 1.8704189044736084, 1.240568989312636], [1.6494649440008793, 2.235214894314338, -0.8029305726284828], [1.9974224573334813, -0.6906780683055276, 1.0535722235492986], [1.9277183224308947, -0.6529932317206242, -1.110224125209533], [-0.3501493572534769, 0.7132786644586329, 1.0706086973199325], [-0.42668432219276026, -1.2491488329668505, 0.15362238828850008], [1.882085252168938, -3.1567731742261307, 0.8497623406176102]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0028', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
