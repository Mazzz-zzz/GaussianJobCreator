import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0052'
logfile = 'conf/5009017845242299296281_0052.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863839, 0.6217394783082123, -1.2501828803165018], [-0.39761971585595596, -0.07566485901595937, -2.6339101198206674], [1.087941409756326, -0.08367235852192238, -3.1385123060458486], [1.3412747433710004, -1.226837997824201, -4.1818625359304695], [0.31359443784925833, -1.3166825080227524, -5.007878938132581], [2.4483457429192925, -0.9821601217258167, -4.877151381993846], [1.5638064896882018, -2.893937332736271, -3.349442713646216], [1.3877310311824553, -3.908088204968372, -4.330730871667522], [2.702177123588236, -2.787608558033008, -2.5138810682889527], [0.2890143247871003, -2.8586915619650206, -2.4252104366949583], [1.350706773188925, 1.0832029916023957, -3.7201890492274705], [1.9101736356828996, -0.26833351767859953, -2.1136883092367644], [-0.7843720924895173, -1.34454935469888, -2.4752210448161533], [-1.1581983763610546, 0.49302227342428917, -3.5614953002063574], [-0.2560445759534288, 1.8735740976390411, -1.3916701657561943], [-2.007602477187449, 0.6427130616946715, -1.0543092166280625], [1.5770424436171646, 0.0, 0.0], [2.292718146893913, 1.3915527243580548, 0.0], [2.341079856722317, 2.05985269289493, 1.4165023767064737], [1.1453916410070326, 2.033172555891593, 1.975951122830792], [3.202217992284641, 1.4422259961598267, 2.1984030628363267], [2.7215555591494947, 3.3152059779287977, 1.2745358845394017], [1.6292323391939727, 2.2122558673105814, -0.8090479336198872], [3.5455868300943787, 1.260039221431076, -0.42809146886198424], [1.9974224573334824, -0.6906780683055259, 1.053572223549299], [1.9277183224308942, -0.6529932317206257, -1.110224125209531], [-0.3501493572534745, 0.7132786644586332, 1.070608697319931], [-0.42668432219275726, -1.2491488329668514, 0.15362238828850358], [0.021434877808575684, -1.9444203672312819, -2.25647452252488]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0052', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
