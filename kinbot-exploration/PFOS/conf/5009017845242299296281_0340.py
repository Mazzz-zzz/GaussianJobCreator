import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0340'
logfile = 'conf/5009017845242299296281_0340.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863849, 0.621739478308216, -1.250182880316501], [-0.397619715855957, -0.075664859015957, -2.6339101198206665], [-1.1233882121466363, 0.4940238789133657, -3.902883790150537], [-2.5854106473305376, -0.05829760723788628, -4.033331227674649], [-3.1611759047412433, -0.09001992508052208, -2.844218110268874], [-3.2994930518117314, 0.7102213573343772, -4.850867135911192], [-2.6135892170537445, -1.7973642989024439, -4.737970384453662], [-3.8849501725471276, -2.369995270099337, -4.45838837594266], [-2.0429277028358066, -1.7259654024627722, -6.032035454465866], [-1.545196823703424, -2.447919188096705, -3.7809653375570327], [-1.1808982350936261, 1.8193192695749123, -3.806570357764953], [-0.4526295328666353, 0.15318954183902195, -4.995744508783356], [0.91867026292174, 0.05926285465343823, -2.81830830472323], [-0.6863393184328426, -1.366668279436514, -2.5232286857846624], [-0.25604457595343155, 1.8735740976390427, -1.3916701657561943], [-2.0076024771874468, 0.6427130616946722, -1.0543092166280594], [1.5770424436171662, 0.0, 0.0], [2.292718146893912, 1.3915527243580592, 0.0], [3.782335574419711, 1.3186147352454682, -0.4807937764724471], [4.419592497958775, 0.34925073248440375, 0.14976318527134808], [3.8473329977166353, 1.1104432852892048, -1.7795740852228796], [4.368010411051574, 2.4684340853352174, -0.20560554420872013], [2.2938967322202157, 1.8704189044736133, 1.2405689893126315], [1.6494649440008762, 2.235214894314341, -0.802930572628488], [1.997422457333486, -0.6906780683055225, 1.0535722235492992], [1.927718322430896, -0.6529932317206203, -1.110224125209533], [-0.35014935725347657, 0.7132786644586356, 1.0706086973199327], [-0.4266843221927533, -1.249148832966848, 0.15362238828850427], [-1.0179922218691675, -3.1018289770861163, -4.260985839582321]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0340', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
