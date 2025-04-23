import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0250'
logfile = 'conf/5009017845242299296281_0250.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863851, -1.393559872884597, 0.08664925740765075], [-0.3976197158559567, -2.2432006455416396, 1.3824827499919141], [-1.1233882121466345, -3.627008449745547, 1.5236046658601572], [-0.3678111583333635, -4.58148330604694, 2.5123747217585426], [0.6688104135818967, -5.129831292172003, 1.903086492896763], [0.05149838568730978, -3.9059088851453563, 3.5785260940644052], [-1.475740760567105, -5.964456490011123, 3.130064368786193], [-2.224617795791914, -6.460477568220277, 2.0275321543384126], [-0.6859569336734569, -6.774825107752932, 3.981475539610524], [-2.428681401794962, -5.102813068078201, 4.041117263817528], [-2.351535285418466, -3.4226112841694123, 1.9914380952175195], [-1.1885301103712396, -4.223582124447545, 0.3402677800844897], [0.9186702629217409, -2.4703580149136957, 1.3578310147309458], [-0.686339318432839, -1.5018460017288933, 2.4451837914307184], [-0.25604457595343316, -2.1420087660532836, -0.9267276815498274], [-2.007602477187449, -1.2344150958913087, -0.02945123045763408], [1.5770424436171657, 0.0, 0.0], [2.292718146893921, 1.3915527243580528, 0.0], [1.60052154700826, 2.440721904563901, -0.9357086002340145], [1.3760692761371707, 1.9145350254105862, -2.1257143081021352], [0.4611709133627646, 2.851946778983259, -0.41882897761344096], [2.4048649076935136, 3.4780800111829073, -1.06893034033068], [3.5400592233304704, 1.2295174124845945, -0.43152105569274457], [2.3088468039522647, 1.896094738758367, 1.2310220414904731], [1.9974224573334851, -0.6906780683055306, 1.0535722235492928], [1.9277183224308914, -0.6529932317206257, -1.110224125209539], [-0.3501493572534819, 0.5705349971623139, -1.1530217920585797], [-0.42668432219275126, 0.7576153073313031, 1.00498342831273], [-2.6831538097969703, -5.603127725524252, 4.829083902582334]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0250', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
