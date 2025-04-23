import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0198'
logfile = 'conf/5009017845242299296281_0198.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863829, 0.621739478308216, -1.2501828803164978], [-0.3976197158559581, -0.07566485901595582, -2.633910119820666], [1.0879414097563231, -0.08367235852191822, -3.1385123060458495], [1.473042089945532, 1.2684473752977068, -3.833533780919367], [1.0004737759019127, 1.2895411745485015, -5.067525751424902], [0.9773040172258539, 2.2971743117836962, -3.1517239144922664], [3.3329235748106685, 1.5027949558458789, -3.920628010705674], [3.609671411865194, 2.509003786807459, -4.886746756507954], [3.806542766823602, 1.5323627643829767, -2.586398878471138], [3.7106516445294835, 0.1003263311188394, -4.529681285724982], [1.894753649284139, -0.2666142022536053, -2.0970816341016434], [1.2598103705638528, -1.0641590715978355, -4.015744696739173], [-0.7843720924895168, -1.344549354698877, -2.475221044816153], [-1.1581983763610573, 0.4930222734242935, -3.5614953002063547], [-0.25604457595342994, 1.8735740976390454, -1.391670165756191], [-2.007602477187446, 0.6427130616946752, -1.054309216628057], [1.5770424436171657, 0.0, 0.0], [2.292718146893915, 1.3915527243580585, 0.0], [2.341079856722319, 2.0598526928949297, 1.4165023767064735], [1.1453916410070362, 2.0331725558915905, 1.9759511228307929], [3.2022179922846434, 1.4422259961598276, 2.1984030628363245], [2.7215555591494955, 3.315205977928798, 1.2745358845394026], [1.6292323391939747, 2.2122558673105797, -0.8090479336198853], [3.5455868300943796, 1.26003922143108, -0.42809146886198307], [1.997422457333486, -0.6906780683055234, 1.0535722235492988], [1.9277183224308956, -0.6529932317206218, -1.1102241252095322], [-0.35014935725347274, 0.713278664458636, 1.0706086973199331], [-0.42668432219275254, -1.249148832966852, 0.15362238828850366], [3.026496535871437, -0.5514585415258024, -4.32214658071228]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0198', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
