import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0173'
logfile = 'conf/5009017845242299296281_0173.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863828, -1.3935598728845986, 0.08664925740765334], [-2.270962283629193, -1.3932545648232175, 0.053676368673278685], [-2.969991788512713, -0.7382819684025457, -1.1888806696804008], [-2.1956407956954176, -1.0418466209202166, -2.5184724117827257], [-2.9912731058630335, -0.8480107816384679, -3.5557278468151385], [-1.1300133907845047, -0.25306365351390975, -2.624111915688928], [-1.584021457030934, -2.814923365908002, -2.5799864107597315], [-1.2617908272916076, -3.127700132838376, -3.9293451929731136], [-0.6920488892418221, -2.9801411848773833, -1.4925836416558833], [-2.9335170727851736, -3.5342731203057456, -2.203236966543831], [-4.2015349880113835, -1.2287045207786709, -1.2974100223711222], [-3.0244357979279366, 0.5781431133077346, -1.0323857637828688], [-2.621205471792927, -2.6819879407832805, 0.09005605494850606], [-2.713377980741127, -0.7960712600074734, 1.1536489463716062], [-0.3710451618282802, -1.949079984012114, 1.2553873021032869], [-0.24552532002049018, -2.1655580448419807, -0.8973245689258783], [1.5770424436171655, 0.0, 0.0], [2.292718146893918, 1.391552724358056, 0.0], [2.34107985672233, 2.059852692894932, 1.4165023767064664], [1.1453916410070473, 2.0331725558915967, 1.975951122830793], [3.2022179922846576, 1.442225996159822, 2.1984030628363245], [2.7215555591495115, 3.3152059779287955, 1.2745358845393973], [1.6292323391939745, 2.212255867310577, -0.8090479336198887], [3.545586830094383, 1.2600392214310727, -0.4280914688619896], [1.997422457333486, -0.6906780683055274, 1.053572223549294], [1.9277183224308931, -0.6529932317206216, -1.1102241252095353], [-0.35014935725347734, 0.5705349971623108, -1.1530217920585788], [-0.4266843221927564, 0.7576153073313028, 1.0049834283127295], [-3.4009948719267906, -3.812129959338628, -3.003375998842822]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0173', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
