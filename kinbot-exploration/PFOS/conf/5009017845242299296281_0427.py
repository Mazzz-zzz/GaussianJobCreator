import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0427'
logfile = 'conf/5009017845242299296281_0427.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863819, -1.3935598728846, 0.08664925740764845], [-2.270962283629193, -1.3932545648232169, 0.05367636867326874], [-3.020318448930584, -0.6573324959322188, 1.2192077454993038], [-4.480058160568045, -0.25598657285097415, 0.8097366435926613], [-4.45376171327607, 0.8374671518423233, 0.06800259968529677], [-5.056700168542038, -1.2407959453510558, 0.1267033098316143], [-5.555811090753356, 0.08585086953294412, 2.3089251331753196], [-5.935868144708009, -1.1611268258091485, 2.877286172364878], [-4.931749541331282, 1.130514850023068, 3.033204859556961], [-6.7987732354635675, 0.7000919263438845, 1.5617332831433688], [-3.086116569243985, -1.4709098648755075, 2.26932870263355], [-2.3639288743924594, 0.44751541508172105, 1.549198080895338], [-2.5961406643712652, -0.7866250332999379, -1.0913507571730194], [-2.6897770261666514, -2.652128594710339, 0.007528686185000967], [-0.3710451618282821, -1.9490799840121147, 1.2553873021032775], [-0.2455253200204908, -2.165558044841978, -0.8973245689258863], [1.5770424436171673, 0.0, 0.0], [2.292718146893918, 1.3915527243580523, 0.0], [3.782335574419717, 1.3186147352454585, -0.4807937764724422], [4.419592497958778, 0.34925073248439364, 0.14976318527135185], [3.847332997716643, 1.1104432852891992, -1.7795740852228734], [4.368010411051581, 2.4684340853351996, -0.20560554420871546], [2.2938967322202215, 1.8704189044735968, 1.2405689893126348], [1.6494649440008855, 2.2352148943143364, -0.8029305726284803], [1.9974224573334836, -0.6906780683055309, 1.0535722235492937], [1.9277183224308951, -0.6529932317206212, -1.1102241252095373], [-0.3501493572534747, 0.5705349971623119, -1.1530217920585784], [-0.42668432219275393, 0.7576153073312978, 1.0049834283127284], [-7.468324495646526, 0.019049892962373038, 1.4073142252166013]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0427', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
