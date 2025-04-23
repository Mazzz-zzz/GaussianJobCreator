import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0021'
logfile = 'conf/5009017845242299296281_0021.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863831, 0.771820394576386, 1.1635336229088467], [-0.39761971585595596, 2.318865504557594, 1.2514273698287413], [1.087941409756327, 2.759867566386749, 1.4967937649483765], [1.4730420899455359, 2.685713952893135, 3.015274540831184], [2.7879690257787098, 2.6561530893205583, 3.144097100354801], [0.9875700214111515, 3.737929380466357, 3.668006955533304], [0.7750378848301916, 1.1500169232449506, 3.8375203032927727], [1.474093054205148, 0.9424904164155907, 5.0584566763046315], [-0.6348786909392101, 1.233472542013578, 3.734918253047024], [1.2673983840267082, 0.08545289608235905, 2.786412696664736], [1.8947536492841404, 1.9494330700686107, 0.8176461448894714], [1.2598103705638566, 4.009816458287681, 1.0862835586981863], [-0.7843720924895167, 2.815878982142089, 0.07319662459686825], [-1.158198376361053, 2.8378342687254494, 2.2077174635201686], [-0.2560445759534292, 0.2684346684142469, 2.318397847306012], [-2.007602477187447, 0.5917020341966387, 1.0837604470856927], [1.5770424436171653, 0.0, 0.0], [2.2927181468939177, 1.391552724358054, 0.0], [2.34107985672233, 2.0598526928949283, 1.416502376706473], [1.1453916410070384, 2.033172555891596, 1.9759511228307893], [3.2022179922846488, 1.4422259961598263, 2.1984030628363254], [2.721555559149512, 3.315205977928792, 1.2745358845394033], [1.6292323391939831, 2.2122558673105748, -0.8090479336198876], [3.5455868300943845, 1.2600392214310667, -0.4280914688619795], [1.9974224573334838, -0.6906780683055265, 1.0535722235493004], [1.9277183224308958, -0.6529932317206288, -1.1102241252095313], [-0.3501493572534749, -1.2838136616209435, 0.08241309473865185], [-0.42668432219275865, 0.4915335256355409, -1.1586058166012296], [1.4019409220005559, 0.5022346949444221, 1.9237693031881153]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0021', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
