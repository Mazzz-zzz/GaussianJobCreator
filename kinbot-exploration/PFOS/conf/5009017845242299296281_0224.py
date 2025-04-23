import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0224'
logfile = 'conf/5009017845242299296281_0224.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863824, 0.771820394576388, 1.1635336229088453], [-0.3466020415139052, 0.2926212946843472, 2.6259377241923785], [-0.7363023803695233, -1.1753096862431975, 3.019500458011147], [0.11232051632788764, -1.6971162941224924, 4.230756914000712], [0.2504583136006516, -0.7374350235871703, 5.128772365440905], [-0.47946595887502164, -2.7474668551256096, 4.79240850290647], [1.8265739031510884, -2.2433886489856847, 3.697291665655142], [1.7145076756829611, -3.5191230914630576, 3.078927856277517], [2.4467690461487464, -1.1216839760476516, 3.0952174318939716], [2.4612831805204625, -2.4290116527713193, 5.126707932351394], [-2.0209682759934693, -1.2016995572775981, 3.3629078892550193], [-0.5312999584009234, -1.9866202042857624, 1.9898504923482956], [0.9809271994195066, 0.4056495998099097, 2.724518672905618], [-0.9145558447621711, 1.1282242700023815, 3.487038975675317], [-2.0119541879597174, 0.6059455318059243, 1.0475948322279907], [-0.4099470658637744, 2.065478369993789, 1.0587059160250698], [1.577042443617164, 0.0, 0.0], [2.292718146893918, 1.3915527243580532, 0.0], [2.341079856722329, 2.059852692894926, 1.4165023767064704], [1.1453916410070413, 2.033172555891602, 1.9759511228307898], [3.202217992284648, 1.4422259961598263, 2.1984030628363254], [2.721555559149514, 3.3152059779287955, 1.2745358845394004], [1.629232339193984, 2.212255867310576, -0.8090479336198888], [3.5455868300943845, 1.2600392214310656, -0.4280914688619817], [1.997422457333482, -0.6906780683055266, 1.0535722235493012], [1.9277183224308951, -0.6529932317206277, -1.110224125209529], [-0.35014935725347834, -1.2838136616209417, 0.08241309473865059], [-0.4266843221927583, 0.4915335256355474, -1.1586058166012296], [1.9952003216854641, -1.8825923575905663, 5.7749016605302135]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0224', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
