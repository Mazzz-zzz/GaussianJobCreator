import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0438'
logfile = 'conf/5009017845242299296281_0438.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863837, 0.7718203945763864, 1.1635336229088453], [-0.3976197158559576, 2.3188655045575954, 1.251427369828738], [-1.123388212146636, 3.132984570832184, 2.37927912429036], [-0.367811158333368, 4.466521985892247, 2.7114935691716617], [0.05889782664840387, 5.030881614919977, 1.5952858849405822], [-1.1736420480321823, 5.307786300766391, 3.3532185650067654], [1.1196427974117333, 4.168630239968735, 3.8162127952032754], [1.9616071617376694, 5.312604247290193, 3.7465673532267636], [0.6301620320434117, 3.6162554763669488, 5.024785288899168], [1.7841101462045712, 3.001651700664524, 2.99329915392575], [-2.351535285418468, 3.4359416226071806, 1.9683492717612094], [-1.1885301103712402, 2.4064716038662994, 3.4875955246991617], [0.9186702629217393, 2.411095160260259, 1.460477289992265], [-0.6863393184328443, 2.8685142811653965, 0.07804489435392932], [-0.2560445759534295, 0.26843466841425334, 2.318397847306014], [-2.007602477187447, 0.5917020341966347, 1.0837604470856923], [1.5770424436171657, 0.0, 0.0], [2.292718146893915, 1.3915527243580539, 0.0], [1.600521547008245, 2.440721904563901, -0.9357086002340216], [1.3760692761371627, 1.9145350254105764, -2.1257143081021397], [0.4611709133627453, 2.8519467789832462, -0.41882897761344984], [2.404864907693492, 3.4780800111829135, -1.0689303403306794], [3.540059223330462, 1.229517412484612, -0.43152105569273913], [2.3088468039522425, 1.8960947387583744, 1.231022041490475], [1.9974224573334836, -0.6906780683055224, 1.0535722235493028], [1.9277183224308976, -0.6529932317206311, -1.1102241252095275], [-0.35014935725347246, -1.2838136616209448, 0.08241309473865537], [-0.42668432219275526, 0.49153352563553443, -1.1586058166012283], [1.1184069582050495, 2.5450746301291125, 2.4600646550000183]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0438', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
