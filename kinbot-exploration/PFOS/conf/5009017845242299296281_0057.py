import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0057'
logfile = 'conf/5009017845242299296281_0057.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863817, -1.393559872884598, 0.0866492574076533], [-0.3466020415139013, -2.420439425248669, -1.0595513872112527], [-0.9873929842445225, -3.849810268767849, -0.9708648395635571], [-0.1847769512407465, -4.7876250073543885, -0.0034475673086363707], [0.8897635883676658, -5.251777005738082, -0.6169844400485595], [0.18128897051813278, -4.122301251114864, 1.0884499008005055], [-1.2119219398679453, -6.254814906871081, 0.5569008084445375], [-0.3308034559149338, -7.240352134917068, 1.0810931368715915], [-2.325919188868468, -5.732942665224553, 1.2583395187151585], [-1.7186561478585216, -6.7307739896164005, -0.8563316649346288], [-2.230506898554333, -3.735942783996664, -0.5118511747699866], [-1.0036314587762272, -4.410925069658703, -2.1730779211988813], [-0.7784534903451347, -1.8458422644301322, -2.1856818928204587], [0.9731582247379494, -2.5497215866483636, -1.121342832061477], [-2.011954187959715, -1.2102165034856978, 0.0009668077395988351], [-0.40994706586377827, -1.949605403411481, 1.2594037813693544], [1.5770424436171682, 0.0, 0.0], [2.292718146893919, 1.3915527243580539, 0.0], [2.3410798567223274, 2.059852692894925, 1.4165023767064744], [1.1453916410070368, 2.033172555891594, 1.9759511228307916], [3.202217992284654, 1.4422259961598174, 2.1984030628363342], [2.7215555591495035, 3.3152059779287897, 1.2745358845394041], [1.6292323391939831, 2.2122558673105726, -0.8090479336198905], [3.5455868300943845, 1.260039221431073, -0.42809146886197874], [1.9974224573334824, -0.6906780683055241, 1.0535722235493028], [1.9277183224308976, -0.6529932317206295, -1.1102241252095306], [-0.35014935725347157, 0.5705349971623058, -1.1530217920585826], [-0.42668432219275454, 0.7576153073313074, 1.0049834283127226], [-1.7082234698932064, -5.993660007613391, -1.4828430562613426]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0057', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
