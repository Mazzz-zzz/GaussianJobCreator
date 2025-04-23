import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0024'
logfile = 'conf/5009017845242299296281_0024.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863849, 0.7718203945763881, 1.1635336229088453], [-0.397619715855957, 2.3188655045575968, 1.2514273698287346], [-1.1233882121466363, 3.1329845708321886, 2.379279124290355], [-2.5854106473305376, 3.522116108662277, 1.9661784049894486], [-3.1611759047412433, 2.5081750999368824, 1.3441495131679213], [-3.2994930518117314, 3.8458634914149803, 3.040503305717411], [-2.613589217053744, 5.001884864766421, 0.8124220495220821], [-3.884950172547126, 5.046075228553262, 0.17671807721631427], [-2.042927702835806, 6.086878641327245, 1.5212878426471084], [-1.545196823703424, 4.498371627201147, -0.22947753452463077], [-1.1808982350936261, 2.3869269963298314, 3.4788618839288867], [-0.4526295328666354, 4.2498468845035005, 2.6305382892183493], [0.9186702629217399, 2.41109516026026, 1.4604772899922642], [-0.6863393184328426, 2.868514281165398, 0.07804489435392553], [-0.25604457595343155, 0.2684346684142553, 2.3183978473060134], [-2.007602477187447, 0.5917020341966408, 1.0837604470856874], [1.5770424436171662, 0.0, 0.0], [2.292718146893917, 1.391552724358052, 0.0], [2.341079856722327, 2.0598526928949275, 1.4165023767064728], [1.1453916410070404, 2.0331725558915945, 1.9759511228307898], [3.202217992284643, 1.4422259961598192, 2.198403062836332], [2.7215555591495106, 3.3152059779287915, 1.274535884539402], [1.6292323391939867, 2.212255867310575, -0.8090479336198894], [3.545586830094387, 1.2600392214310632, -0.4280914688619758], [1.997422457333482, -0.6906780683055244, 1.053572223549307], [1.9277183224308967, -0.6529932317206344, -1.1102241252095268], [-0.35014935725347657, -1.2838136616209441, 0.08241309473865532], [-0.42668432219275265, 0.4915335256355354, -1.1586058166012305], [-1.9826640778384341, 4.076234753571363, -0.9820632427184887]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0024', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
