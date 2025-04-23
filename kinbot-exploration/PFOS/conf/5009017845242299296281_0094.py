import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0094'
logfile = 'conf/5009017845242299296281_0094.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863817, 0.6217394783082225, -1.2501828803164972], [-0.3466020415139013, 2.127818130564332, -1.5663863369811124], [-0.9873929842445225, 2.765698749087088, -2.848601072721357], [-0.18477695124074647, 2.3967981845477584, -4.144481096508225], [0.8897635883676658, 3.1602127016908375, -4.239680081955852], [0.18128897051813278, 1.1185253607175698, -4.114242555918093], [-1.211921939867945, 2.645117205934515, -5.695279009542206], [-0.3308034559149338, 2.6839219470707723, -6.810875449618855], [-2.325919188868468, 1.7767173428191052, -5.594043745881696], [-1.7186561478585216, 4.106991970706649, -5.4008554296719975], [-2.230506898554333, 2.3112475123060685, -2.979495770641263], [-1.0036314587762272, 4.087403218990677, -2.733434203914616], [-0.7784534903451347, 2.8157771759892434, -0.5057053459652403], [0.9731582247379494, 2.245972172241018, -1.647452250584297], [-2.011954187959715, 0.6042709716797878, -1.0485616399675883], [-0.40994706586377827, -0.11587296658229075, -2.318109697394429], [1.5770424436171682, 0.0, 0.0], [2.2927181468939257, 1.3915527243580657, 0.0], [3.782335574419716, 1.3186147352454634, -0.480793776472444], [4.419592497958777, 0.3492507324843893, 0.14976318527134663], [3.847332997716635, 1.110443285289206, -1.779574085222875], [4.368010411051582, 2.4684340853352014, -0.20560554420871197], [2.2938967322202206, 1.8704189044735995, 1.2405689893126397], [1.6494649440008835, 2.2352148943143417, -0.8029305726284779], [1.9974224573334887, -0.6906780683055331, 1.053572223549291], [1.927718322430894, -0.6529932317206205, -1.110224125209535], [-0.3501493572534714, 0.7132786644586308, 1.0706086973199378], [-0.42668432219275454, -1.2491488329668499, 0.15362238828849772], [-1.7082234698932062, 4.281009760354404, -4.44924030010933]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0094', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
