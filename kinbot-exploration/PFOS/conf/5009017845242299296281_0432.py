import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0432'
logfile = 'conf/5009017845242299296281_0432.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863837, -1.3935598728845955, 0.08664925740765314], [-0.397619715855957, -2.2432006455416342, 1.382482749991919], [1.0879414097563251, -2.676195207864826, 1.6417185410974655], [1.3412747433709995, -3.0081801923381053, 3.1534041404090263], [2.433578518280429, -3.742064612588835, 3.275010324909891], [1.4798103010878627, -1.8878195086267913, 3.856610993472195], [-0.08698699647460374, -3.9648101239617124, 3.9061264713717923], [-0.4971995188277743, -4.9579059002113635, 2.9745055862958147], [0.2629292866236743, -4.235299785711872, 5.251396224751005], [-1.155974242147192, -2.808206044490851, 3.914838632395782], [1.3507067731889257, -3.7633797193128675, 0.9220132164307515], [1.9101736356828996, -1.6963410126419187, 1.2892277976148874], [-0.7843720924895171, -1.4713296274432097, 2.4020244202192793], [-1.1581983763610553, -3.330856542149738, 1.35377783668618], [-0.2560445759534292, -2.1420087660532836, -0.9267276815498219], [-2.0076024771874477, -1.234415095891308, -0.029451230457632965], [1.577042443617163, 0.0, 0.0], [2.2927181468939146, 1.391552724358058, 0.0], [3.782335574419718, 1.3186147352454574, -0.4807937764724439], [4.419592497958776, 0.34925073248439253, 0.1497631852713508], [3.847332997716636, 1.1104432852891986, -1.7795740852228779], [4.368010411051581, 2.468434085335205, -0.2056055442087159], [2.2938967322202193, 1.8704189044736033, 1.2405689893126353], [1.6494649440008824, 2.2352148943143386, -0.802930572628484], [1.997422457333483, -0.6906780683055257, 1.0535722235492988], [1.9277183224308931, -0.6529932317206223, -1.1102241252095302], [-0.35014935725347457, 0.5705349971623117, -1.153021792058578], [-0.4266843221927585, 0.7576153073313049, 1.0049834283127286], [-0.7178336349735206, -1.945713562296458, 3.9249338391614454]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0432', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
