import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0493'
logfile = 'conf/5009017845242299296281_0493.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863842, -1.3935598728845981, 0.08664925740765331], [-0.34660204151390284, -2.4204394252486665, -1.0595513872112554], [1.1624535490467065, -2.799902192493113, -1.2603960465750377], [1.924360652536601, -1.723803915021475, -2.1095693044418975], [1.6581688080331816, -1.8936103403721951, -3.3928501088282457], [1.5600180827406167, -0.4992730224820855, -1.739868580065883], [3.7832317877839663, -1.8440144326786716, -1.8819493028959793], [4.116829961099136, -1.257966558680993, -0.6298745623339073], [4.161464969515624, -3.1518513587361334, -2.2716402146782992], [4.217392380657199, -0.8594959572251144, -3.0320669407468097], [1.7462981079693096, -2.8864928243592156, -0.06849875821689183], [1.2535978601361863, -3.964626716020741, -1.8893151990408852], [-1.011797100425898, -3.534969835652361, -0.7430628651027207], [-0.8215282967123789, -1.9466956411138827, -2.2051128650057947], [-2.0119541879597187, -1.2102165034856953, 0.0009668077395945335], [-0.40994706586378143, -1.9496054034114758, 1.2594037813693562], [1.577042443617166, 0.0, 0.0], [2.292718146893918, 1.3915527243580557, 0.0], [2.34107985672232, 2.059852692894929, 1.4165023767064728], [1.1453916410070355, 2.033172555891595, 1.9759511228307924], [3.202217992284644, 1.4422259961598267, 2.198403062836328], [2.7215555591494973, 3.3152059779287946, 1.2745358845394008], [1.6292323391939758, 2.212255867310575, -0.8090479336198898], [3.5455868300943822, 1.2600392214310707, -0.4280914688619808], [1.9974224573334802, -0.6906780683055239, 1.0535722235493024], [1.9277183224308965, -0.6529932317206275, -1.1102241252095284], [-0.3501493572534737, 0.5705349971623084, -1.1530217920585826], [-0.4266843221927596, 0.7576153073313069, 1.004983428312722], [5.034132539369231, -1.1685200040623127, -3.448482644484556]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0493', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
