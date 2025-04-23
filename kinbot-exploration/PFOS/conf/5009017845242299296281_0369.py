import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0369'
logfile = 'conf/5009017845242299296281_0369.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863858, 0.7718203945763842, 1.1635336229088462], [-2.2709622836291965, 0.7431123812655637, 1.1797556627388934], [-3.020318448930585, 1.384531128059256, -0.04033723253933168], [-3.1314232385001857, 2.942629924600787, 0.09955209854967853], [-1.9982160731709881, 3.433107501124375, 0.5704822530854748], [-3.396984742281267, 3.495707533466796, -1.080482768237572], [-4.510308809776061, 3.441050521602605, 1.270839774966773], [-5.748558712762754, 3.2997461390756326, 0.5858571058122058], [-4.223125717516355, 2.843275764194895, 2.522058374003057], [-4.168531436650326, 4.974230245459327, 1.3852547712867498], [-2.345438986660112, 1.1086318389066705, -1.1527260342434005], [-4.247573342377324, 0.8876550485557878, -0.12629879149041923], [-2.5961406643712657, -0.5518249635012427, 1.226912640677034], [-2.689777026166658, 1.3325843308485121, 2.2930463940297616], [-0.37104516182828684, 2.061737287215913, 1.0602591291106052], [-0.24552532002049313, 0.3056731502912673, 2.3240905646658523], [1.5770424436171646, 0.0, 0.0], [2.2927181468939137, 1.3915527243580552, 0.0], [3.7823355744197142, 1.3186147352454636, -0.4807937764724375], [4.419592497958774, 0.34925073248439875, 0.1497631852713578], [3.8473329977166415, 1.1104432852891986, -1.7795740852228714], [4.368010411051579, 2.468434085335212, -0.20560554420871568], [2.2938967322202113, 1.8704189044736121, 1.2405689893126344], [1.6494649440008788, 2.2352148943143377, -0.8029305726284879], [1.997422457333483, -0.6906780683055247, 1.0535722235493035], [1.9277183224308976, -0.6529932317206257, -1.1102241252095246], [-0.3501493572534754, -1.2838136616209455, 0.08241309473865198], [-0.4266843221927567, 0.4915335256355403, -1.1586058166012294], [-3.2331902639829475, 5.125529899253409, 1.1897868911350307]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0369', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
