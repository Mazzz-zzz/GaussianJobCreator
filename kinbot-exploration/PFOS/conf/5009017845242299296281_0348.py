import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0348'
logfile = 'conf/5009017845242299296281_0348.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863863, 0.7718203945763877, 1.1635336229088473], [-0.3976197158559574, 2.318865504557596, 1.251427369828742], [1.0879414097563262, 2.7598675663867502, 1.4967937649483756], [1.473042089945535, 2.685713952893136, 3.0152745408311827], [2.7879690257787106, 2.6561530893205587, 3.1440971003548013], [0.9875700214111528, 3.7379293804663556, 3.6680069555333064], [0.7750378848301898, 1.1500169232449544, 3.8375203032927727], [-0.6091697343464534, 1.3663644674714959, 4.081639157399683], [1.2916422504512648, 0.03553178997703076, 3.132963996959468], [1.5508075901178944, 1.2513017732964247, 5.204438506155287], [1.8947536492841404, 1.9494330700686107, 0.8176461448894714], [1.2598103705638566, 4.009816458287681, 1.0862835586981863], [-0.7843720924895178, 2.815878982142091, 0.0731966245968683], [-1.1581983763610542, 2.8378342687254516, 2.2077174635201664], [-0.2560445759534314, 0.2684346684142489, 2.318397847306012], [-2.0076024771874486, 0.591702034196641, 1.0837604470856925], [1.5770424436171633, 0.0, 0.0], [2.2927181468939164, 1.391552724358054, 0.0], [2.341079856722329, 2.0598526928949283, 1.4165023767064726], [1.145391641007037, 2.0331725558915967, 1.9759511228307902], [3.2022179922846483, 1.4422259961598256, 2.198403062836325], [2.7215555591495115, 3.3152059779287923, 1.2745358845394026], [1.6292323391939822, 2.212255867310575, -0.8090479336198877], [3.5455868300943836, 1.2600392214310658, -0.4280914688619803], [1.9974224573334822, -0.6906780683055262, 1.0535722235493008], [1.9277183224308945, -0.6529932317206295, -1.11022412520953], [-0.35014935725347546, -1.2838136616209435, 0.08241309473865185], [-0.426684322192762, 0.4915335256355432, -1.158605816601229], [1.7803408339398903, 0.368486163200742, 5.526794637240062]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0348', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
