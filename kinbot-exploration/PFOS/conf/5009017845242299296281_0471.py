import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0471'
logfile = 'conf/5009017845242299296281_0471.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863824, 0.7718203945763835, 1.163533622908852], [-0.3976197158559557, 2.3188655045575906, 1.2514273698287526], [-0.7491833517666568, 3.196572017163506, -0.0006164121359253573], [0.3751315092500395, 3.1292034692660815, -1.0919288014421222], [0.8222268376054805, 1.8902897582872573, -1.200097776607916], [-0.09608992658046996, 3.5343107186271188, -2.26782782579837], [1.8294372022557681, 4.235362410296145, -0.6640612055840428], [2.112852701691653, 4.0815510595152755, 0.7210195613246831], [2.788205081300309, 4.0704230018413, -1.6930976079037108], [1.128311689997794, 5.62834272973815, -0.8844505780983333], [-1.8820373528306351, 2.7491646570443438, -0.5349009981858436], [-0.8963077255402097, 4.463477707111217, 0.3650804471736064], [-1.1305489784219276, 2.7568196252985233, 2.278967330580903], [0.887865254859781, 2.484225820671751, 1.5386269105814956], [-0.25604457595342717, 0.26843466841424124, 2.318397847306018], [-2.007602477187446, 0.5917020341966374, 1.0837604470856996], [1.5770424436171644, 0.0, 0.0], [2.292718146893916, 1.3915527243580539, 0.0], [2.3410798567223283, 2.05985269289493, 1.416502376706471], [1.1453916410070395, 2.033172555891594, 1.975951122830792], [3.2022179922846488, 1.4422259961598267, 2.198403062836321], [2.721555559149507, 3.315205977928799, 1.274535884539398], [1.6292323391939743, 2.2122558673105774, -0.8090479336198897], [3.545586830094382, 1.260039221431074, -0.4280914688619871], [1.9974224573334838, -0.6906780683055265, 1.0535722235493028], [1.9277183224308951, -0.6529932317206233, -1.1102241252095308], [-0.3501493572534753, -1.2838136616209441, 0.08241309473865305], [-0.42668432219275954, 0.4915335256355467, -1.1586058166012259], [0.38543370612337463, 5.537403989715657, -1.4974961174620889]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0471', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
