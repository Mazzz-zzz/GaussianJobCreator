import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0260'
logfile = 'conf/5009017845242299296281_0260.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863852, -1.3935598728845988, 0.08664925740765086], [-0.3466020415139046, -2.420439425248664, -1.0595513872112607], [1.1624535490467056, -2.7999021924931085, -1.2603960465750463], [1.322955526797161, -4.1716963707831285, -2.003366440702432], [2.5379668391891146, -4.265983898549034, -2.514599573377115], [1.1168882370515751, -5.1808849811579964, -1.1619541390693058], [0.08826989713527927, -4.354401773287279, -3.4047348658202665], [-1.1777533685614712, -4.676778155818564, -2.842872278962436], [0.31035495620094256, -3.273264341132757, -4.291964740758064], [0.6897122389832802, -5.64883116479621, -4.07051432603154], [1.7500514572469157, -1.8550381832545335, -1.9890607533025928], [1.7633282149005813, -2.89451491443956, -0.08123199086868528], [-1.0117971004258985, -3.53496983565236, -0.7430628651027305], [-0.821528296712383, -1.9466956411138774, -2.2051128650057987], [-2.01195418795972, -1.2102165034856998, 0.0009668077395977617], [-0.40994706586377916, -1.9496054034114847, 1.259403781369352], [1.5770424436171668, 0.0, 0.0], [2.2927181468939173, 1.391552724358056, 0.0], [3.7823355744197174, 1.3186147352454587, -0.4807937764724436], [4.419592497958776, 0.3492507324843883, 0.1497631852713539], [3.84733299771664, 1.1104432852891977, -1.7795740852228756], [4.368010411051584, 2.468434085335198, -0.20560554420871502], [2.293896732220221, 1.8704189044736046, 1.240568989312633], [1.6494649440008846, 2.235214894314339, -0.8029305726284793], [1.9974224573334824, -0.6906780683055276, 1.0535722235492981], [1.927718322430894, -0.6529932317206246, -1.1102241252095306], [-0.3501493572534747, 0.5705349971623094, -1.153021792058581], [-0.42668432219275654, 0.7576153073313029, 1.0049834283127295], [0.592277074907428, -5.608363491177203, -5.032195740476766]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0260', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
