import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0126'
logfile = 'conf/5009017845242299296281_0126.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863841, 0.7718203945763837, 1.163533622908848], [-0.39761971585595707, 2.3188655045575914, 1.2514273698287461], [1.0879414097563251, 2.759867566386746, 1.4967937649483813], [2.100632021959493, 1.8490216187922435, 0.7193604606011571], [3.2636602341713647, 2.466189412327759, 0.6054691295054131], [2.2719309703990027, 0.695851958252546, 1.3594826714541275], [1.48866165575111, 1.4416571655172146, -1.007257564384598], [0.4899953317481413, 0.43496974095033736, -0.8991973481563169], [1.3031116920812955, 2.6743055280880643, -1.6791995238435466], [2.80710863078977, 0.7870438926650963, -1.5673661849451563], [1.2395907752722293, 4.0125250124967184, 1.0762663043113119], [1.3717373297746493, 2.6844668573495607, 2.7907042758982867], [-0.7843720924895173, 2.815878982142089, 0.07319662459687326], [-1.1581983763610537, 2.8378342687254468, 2.207717463520175], [-0.25604457595342883, 0.26843466841424085, 2.3183978473060143], [-2.007602477187449, 0.5917020341966351, 1.0837604470856932], [1.5770424436171655, 0.0, 0.0], [2.2927181468939146, 1.3915527243580554, 0.0], [2.341079856722321, 2.059852692894933, 1.416502376706474], [1.1453916410070322, 2.033172555891592, 1.9759511228307922], [3.202217992284639, 1.442225996159832, 2.198403062836326], [2.721555559149496, 3.315205977928798, 1.274535884539405], [1.6292323391939703, 2.2122558673105757, -0.8090479336198851], [3.545586830094379, 1.2600392214310787, -0.4280914688619816], [1.9974224573334836, -0.6906780683055276, 1.053572223549301], [1.9277183224308962, -0.6529932317206264, -1.1102241252095315], [-0.3501493572534735, -1.2838136616209448, 0.08241309473864851], [-0.4266843221927592, 0.49153352563554376, -1.1586058166012283], [2.769835737892855, -0.17530045946234793, -1.4753527833152342]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0126', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
