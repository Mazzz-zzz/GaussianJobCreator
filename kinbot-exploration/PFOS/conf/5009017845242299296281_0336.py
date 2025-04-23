import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0336'
logfile = 'conf/5009017845242299296281_0336.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863845, 0.7718203945763843, 1.1635336229088487], [-0.397619715855956, 2.3188655045575914, 1.2514273698287477], [-0.7491833517666574, 3.1965720171635055, -0.000616412135929616], [-2.081987459715323, 2.729224887216827, -0.6823825222996486], [-2.976543962833394, 2.4254784423769804, 0.24174339068695233], [-2.563171468208477, 3.690357701799604, -1.4656996376850224], [-1.821073350728775, 1.2133943698125083, -1.757484843640507], [-3.089985411108095, 0.6224049358188903, -2.0085320100092967], [-0.9027569326874728, 1.5811305269912916, -2.7708070267131832], [-1.047509876275078, 0.3240325775280425, -0.7128902891722709], [-0.8953137364216059, 4.460279692037433, 0.38755928990030625], [0.22590414679516443, 3.113269857861802, -0.8965488436803655], [-1.1305489784219298, 2.756819625298524, 2.2789673305808975], [0.8878652548597805, 2.4842258206717514, 1.5386269105814916], [-0.2560445759534273, 0.26843466841424285, 2.318397847306014], [-2.0076024771874477, 0.5917020341966369, 1.083760447085694], [1.5770424436171648, 0.0, 0.0], [2.2927181468939155, 1.3915527243580552, 0.0], [3.782335574419714, 1.318614735245464, -0.4807937764724448], [4.4195924979587735, 0.34925073248439187, 0.1497631852713477], [3.8473329977166357, 1.1104432852892003, -1.779574085222876], [4.3680104110515785, 2.468434085335205, -0.20560554420872013], [2.2938967322202157, 1.8704189044736044, 1.2405689893126348], [1.6494649440008804, 2.235214894314339, -0.802930572628489], [1.9974224573334813, -0.6906780683055285, 1.0535722235492995], [1.9277183224308942, -0.6529932317206258, -1.1102241252095337], [-0.35014935725347784, -1.2838136616209448, 0.08241309473864851], [-0.426684322192759, 0.4915335256355461, -1.1586058166012265], [-0.37293782431117783, -0.2075023335756361, -1.1582952416404948]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0336', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
