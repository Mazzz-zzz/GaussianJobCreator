import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0114'
logfile = 'conf/5009017845242299296281_0114.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863835, 0.7718203945763854, 1.163533622908846], [-2.270962283629194, 0.743112381265571, 1.1797556627388939], [-2.9699917885127154, -0.660459877810206, 1.2338112746327814], [-3.0567398040787155, -1.3216748454235774, -0.18564450396738677], [-4.057927071959131, -0.792630756448847, -0.866994890796401], [-1.9208515428855029, -1.1373859158669188, -0.8524994275297276], [-3.34387825427727, -3.1729541587223333, -0.07619247938739071], [-4.321509485593893, -3.4127867354508172, 0.9283713655244802], [-3.4371870557562354, -3.653830249844668, -1.4048792085306752], [-1.9318779018872667, -3.5894681757511417, 0.4835518446374543], [-2.265942241823496, -1.4610539209006583, 2.029018773642182], [-4.202897212097877, -0.5310943056708307, 1.7067156695461143], [-2.6212054717929294, 1.4189848017416575, 2.277641661887574], [-2.713377980741128, 1.3971249246107, 0.1125934612033489], [-0.37104516182828085, 2.0617372872159145, 1.0602591291106072], [-0.2455253200204934, 0.3056731502912649, 2.32409056466585], [1.5770424436171644, 0.0, 0.0], [2.292718146893919, 1.3915527243580508, 0.0], [3.7823355744197174, 1.3186147352454514, -0.4807937764724357], [4.419592497958776, 0.349250732484386, 0.14976318527135246], [3.847332997716645, 1.1104432852891912, -1.7795740852228679], [4.368010411051584, 2.4684340853351987, -0.20560554420871258], [2.293896732220217, 1.870418904473603, 1.2405689893126337], [1.6494649440008837, 2.235214894314337, -0.8029305726284883], [1.99742245733348, -0.6906780683055298, 1.053572223549299], [1.9277183224308962, -0.6529932317206304, -1.1102241252095284], [-0.35014935725347784, -1.2838136616209423, 0.08241309473865063], [-0.42668432219275365, 0.4915335256355465, -1.1586058166012292], [-1.6621905944176811, -4.442079588468721, 0.11438396712163197]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0114', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
