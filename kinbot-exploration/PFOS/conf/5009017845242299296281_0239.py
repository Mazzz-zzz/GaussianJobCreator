import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0239'
logfile = 'conf/5009017845242299296281_0239.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863828, -1.3935598728845953, 0.08664925740765188], [-0.397619715855957, -2.2432006455416342, 1.382482749991919], [1.0879414097563236, -2.6761952078648275, 1.6417185410974664], [1.4730420899455332, -3.9541613281908345, 0.8182592400881731], [1.000473775901915, -5.033376622340043, 1.4169874593274274], [0.9773040172258548, -3.8780601315570844, -0.41354935367958306], [3.332923574810669, -4.146760933982902, 0.6588553969111904], [3.79171989599616, -3.232379838169373, -0.3289928480614424], [3.852849871526595, -4.220679515015701, 1.973986364157174], [3.3600157795274916, -5.6003105579497925, 0.053048387279991324], [1.8947536492841397, -1.6828188678150033, 1.2794354892121673], [1.2598103705638541, -2.945657386689844, 2.9294611380409767], [-0.7843720924895167, -1.4713296274432102, 2.4020244202192798], [-1.158198376361056, -3.330856542149738, 1.35377783668618], [-0.2560445759534292, -2.142008766053283, -0.9267276815498228], [-2.007602477187446, -1.2344150958913058, -0.029451230457632913], [1.5770424436171657, 0.0, 0.0], [2.2927181468939186, 1.3915527243580552, 0.0], [3.7823355744197182, 1.3186147352454591, -0.48079377647244165], [4.419592497958777, 0.34925073248439353, 0.14976318527135585], [3.8473329977166415, 1.1104432852891986, -1.7795740852228707], [4.368010411051583, 2.468434085335202, -0.2056055442087179], [2.293896732220222, 1.8704189044736026, 1.2405689893126306], [1.649464944000885, 2.2352148943143373, -0.8029305726284838], [1.9974224573334836, -0.6906780683055255, 1.0535722235493017], [1.9277183224308962, -0.652993231720627, -1.1102241252095344], [-0.35014935725347246, 0.5705349971623118, -1.1530217920585781], [-0.42668432219275676, 0.7576153073313046, 1.0049834283127281], [3.41899316549829, -5.561210850041774, -0.9118118448166428]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0239', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
