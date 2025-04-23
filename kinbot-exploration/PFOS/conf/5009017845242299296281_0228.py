import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0228'
logfile = 'conf/5009017845242299296281_0228.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863814, 0.7718203945763792, 1.1635336229088522], [-2.270962283629191, 0.7431123812655601, 1.1797556627389065], [-3.020318448930585, 1.384531128059258, -0.040337232539313886], [-3.131423238500184, 2.94262992460079, 0.09955209854970477], [-4.113036047167058, 3.2500928357049923, 0.9292537287960602], [-1.9898347010198012, 3.4490347131101102, 0.5569685847209446], [-3.4853445173942257, 3.7672726941604275, -1.5485884786080866], [-4.470675829029668, 2.9990263074796455, -2.2278164988433344], [-3.596447598222816, 5.156554403802908, -1.2977713971896312], [-2.089155928253522, 3.516478541160051, -2.2329146772758834], [-2.3454389866601137, 1.1086318389066785, -1.1527260342433856], [-4.247573342377323, 0.8876550485557944, -0.12629879149040085], [-2.596140664371264, -0.5518249635012475, 1.226912640677039], [-2.689777026166647, 1.332584330848503, 2.2930463940297785], [-0.3710451618282818, 2.0617372872159083, 1.0602591291106176], [-0.24552532002048805, 0.305673150291253, 2.3240905646658554], [1.5770424436171642, 0.0, 0.0], [2.292718146893915, 1.3915527243580539, 0.0], [2.341079856722324, 2.0598526928949292, 1.416502376706471], [1.1453916410070413, 2.0331725558915923, 1.9759511228307935], [3.2022179922846488, 1.4422259961598298, 2.1984030628363227], [2.721555559149502, 3.3152059779288, 1.2745358845393988], [1.6292323391939711, 2.212255867310577, -0.8090479336198895], [3.545586830094379, 1.260039221431076, -0.42809146886198846], [1.997422457333485, -0.6906780683055309, 1.0535722235492933], [1.9277183224308942, -0.6529932317206217, -1.1102241252095344], [-0.35014935725347673, -1.2838136616209463, 0.08241309473864745], [-0.42668432219275976, 0.4915335256355476, -1.1586058166012252], [-1.8506310768882523, 4.267061492004583, -2.794780778654752]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0228', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
