import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0445'
logfile = 'conf/5009017845242299296281_0445.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863818, 0.771820394576383, 1.1635336229088467], [-0.39761971585595707, 2.31886550455759, 1.2514273698287484], [-0.7491833517666582, 3.196572017163505, -0.0006164121359296159], [-0.9314389667275274, 4.706813247194561, 0.3808799534188684], [-0.8151262318136054, 5.459876298882798, -0.6988743385178852], [-2.125609869618872, 4.9042934965067975, 0.9321126142547967], [0.35577392705740385, 5.272324182776366, 1.623837587019553], [-0.019766177887423766, 4.78499014641773, 2.905986632331535], [1.6233182718939367, 5.052902781653188, 1.0317714131419502], [0.06381288580124321, 6.8189191742867665, 1.5659989110961268], [0.24079307197239344, 3.10180301143506, -0.8837681932585597], [-1.879565640212353, 2.7679819246717683, -0.5474035883527842], [-1.1305489784219334, 2.756819625298523, 2.278967330580899], [0.8878652548597799, 2.4842258206717514, 1.5386269105814927], [-0.25604457595343105, 0.26843466841424096, 2.318397847306015], [-2.007602477187446, 0.5917020341966339, 1.0837604470856954], [1.5770424436171646, 0.0, 0.0], [2.2927181468939155, 1.3915527243580559, 0.0], [1.600521547008257, 2.440721904563903, -0.9357086002340288], [1.3760692761371667, 1.9145350254105797, -2.1257143081021432], [0.4611709133627564, 2.8519467789832547, -0.4188289776134497], [2.404864907693501, 3.4780800111829135, -1.068930340330683], [3.5400592233304655, 1.2295174124846058, -0.43152105569274735], [2.3088468039522527, 1.8960947387583713, 1.2310220414904696], [1.9974224573334818, -0.6906780683055265, 1.0535722235493], [1.9277183224308945, -0.6529932317206251, -1.1102241252095326], [-0.35014935725347535, -1.2838136616209446, 0.08241309473865077], [-0.4266843221927569, 0.491533525635545, -1.158605816601226], [0.8858755163704096, 7.31991360227724, 1.661811234104188]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0445', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
