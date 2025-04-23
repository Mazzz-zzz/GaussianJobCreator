import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0278'
logfile = 'conf/5009017845242299296281_0278.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863808, 0.6217394783082213, -1.2501828803164974], [-0.3466020415139, 2.1278181305643318, -1.5663863369811133], [-0.9873929842445205, 2.7656987490870835, -2.8486010727213604], [-0.1847769512407454, 2.3967981845477544, -4.144481096508227], [-0.9430454328354119, 2.5834030116378752, -5.210646886338207], [0.908522328044296, 3.1479156941145856, -4.24153372805386], [0.3705079163707562, 0.604253475955705, -4.1342090657780455], [-0.7088075214977008, -0.1930227891592658, -3.663189161237765], [1.0503640582700182, 0.3774067173772625, -5.355557888224902], [1.455146591050225, 0.7086053523466976, -2.9970141567618946], [-2.230506898554333, 2.3112475123060667, -2.9794957706412637], [-1.0036314587762256, 4.087403218990675, -2.7334342039146184], [-0.7784534903451326, 2.8157771759892434, -0.5057053459652429], [0.9731582247379469, 2.245972172241019, -1.6474522505842975], [-2.0119541879597156, 0.6042709716797867, -1.0485616399675894], [-0.40994706586377827, -0.1158729665822928, -2.3181096973944286], [1.577042443617169, 0.0, 0.0], [2.292718146893921, 1.391552724358053, 0.0], [1.6005215470082552, 2.4407219045638993, -0.9357086002340261], [1.3760692761371676, 1.9145350254105789, -2.1257143081021397], [0.4611709133627657, 2.851946778983255, -0.41882897761344096], [2.404864907693506, 3.4780800111829073, -1.068930340330686], [3.5400592233304717, 1.229517412484597, -0.43152105569275034], [2.308846803952263, 1.8960947387583684, 1.2310220414904693], [1.9974224573334878, -0.6906780683055352, 1.0535722235492915], [1.9277183224308934, -0.6529932317206231, -1.1102241252095353], [-0.3501493572534682, 0.7132786644586321, 1.0706086973199376], [-0.4266843221927556, -1.24914883296685, 0.15362238828849664], [2.201414134241506, 0.12386882327878275, -3.18970648211207]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0278', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
