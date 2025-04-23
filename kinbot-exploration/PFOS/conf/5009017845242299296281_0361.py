import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0361'
logfile = 'conf/5009017845242299296281_0361.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863849, 0.7718203945763893, 1.1635336229088449], [-0.39761971585595673, 2.3188655045575977, 1.2514273698287377], [-1.1233882121466356, 3.132984570832186, 2.3792791242903593], [-1.2235971086413773, 2.315167000798345, 3.713714494853222], [-1.448535500070832, 3.1339929795330352, 4.726349019095459], [-2.2078139831012678, 1.4238504083556371, 3.6377413639598877], [0.3571901289937569, 1.3696231845390576, 4.072492883247369], [0.3489316147660549, 0.993465718949593, 5.443940427207419], [0.5452803995777111, 0.4618532063764237, 3.0020105207087497], [1.378733659895341, 2.55272248798352, 3.8794450071353697], [-0.434359585305253, 4.245343424397676, 2.617793033305205], [-2.353319567400157, 3.442949995028898, 1.989894868502084], [0.9186702629217396, 2.4110951602602593, 1.4604772899922667], [-0.686339318432842, 2.8685142811654005, 0.07804489435392942], [-0.25604457595343166, 0.2684346684142534, 2.3183978473060143], [-2.0076024771874486, 0.5917020341966399, 1.0837604470856903], [1.5770424436171633, 0.0, 0.0], [2.292718146893917, 1.391552724358056, 0.0], [3.7823355744197156, 1.3186147352454602, -0.48079377647244526], [4.419592497958774, 0.34925073248439287, 0.14976318527135313], [3.8473329977166375, 1.1104432852891877, -1.7795740852228783], [4.3680104110515785, 2.468434085335206, -0.20560554420872468], [2.2938967322202157, 1.8704189044736124, 1.2405689893126264], [1.6494649440008824, 2.2352148943143333, -0.8029305726284943], [1.997422457333482, -0.6906780683055234, 1.0535722235493061], [1.9277183224308971, -0.6529932317206308, -1.1102241252095262], [-0.3501493572534757, -1.2838136616209443, 0.08241309473865534], [-0.426684322192754, 0.49153352563553815, -1.15860581660123], [1.589194579797825, 2.957571159705877, 4.732537112886234]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0361', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
