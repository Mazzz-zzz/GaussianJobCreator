import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0088'
logfile = 'conf/5009017845242299296281_0088.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863828, -1.3935598728845986, 0.08664925740765334], [-2.270962283629193, -1.3932545648232175, 0.053676368673278685], [-2.969991788512713, -0.7382819684025457, -1.1888806696804008], [-2.1956407956954176, -1.0418466209202166, -2.5184724117827257], [-2.9912731058630335, -0.8480107816384679, -3.5557278468151385], [-1.1300133907845047, -0.25306365351390975, -2.624111915688928], [-1.584021457030934, -2.814923365908002, -2.5799864107597315], [-1.2617908272916074, -3.1277001328383753, -3.9293451929731127], [-0.6920488892418222, -2.9801411848773824, -1.4925836416558846], [-2.9335170727851736, -3.5342731203057456, -2.203236966543831], [-4.2015349880113835, -1.2287045207786709, -1.2974100223711222], [-3.0244357979279366, 0.5781431133077346, -1.0323857637828688], [-2.621205471792927, -2.6819879407832805, 0.09005605494850606], [-2.713377980741127, -0.7960712600074734, 1.1536489463716062], [-0.3710451618282802, -1.949079984012114, 1.2553873021032869], [-0.24552532002049018, -2.1655580448419807, -0.8973245689258783], [1.5770424436171655, 0.0, 0.0], [2.292718146893918, 1.391552724358056, 0.0], [3.7823355744197116, 1.318614735245467, -0.48079377647244403], [4.419592497958776, 0.34925073248439653, 0.14976318527134436], [3.8473329977166313, 1.110443285289201, -1.7795740852228863], [4.368010411051576, 2.4684340853352094, -0.20560554420872051], [2.2938967322202184, 1.870418904473604, 1.2405689893126366], [1.6494649440008777, 2.2352148943143404, -0.8029305726284848], [1.9974224573334858, -0.6906780683055275, 1.0535722235492946], [1.9277183224308925, -0.6529932317206215, -1.1102241252095357], [-0.35014935725347734, 0.5705349971623108, -1.1530217920585788], [-0.4266843221927564, 0.7576153073313028, 1.0049834283127295], [-3.505752722677643, -2.9349998746787502, -1.7038570262540387]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0088', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
