import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0206'
logfile = 'conf/5009017845242299296281_0206.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863845, 0.6217394783082119, -1.2501828803165036], [-0.3976197158559556, -0.07566485901596118, -2.633910119820669], [1.087941409756327, -0.08367235852192448, -3.1385123060458486], [2.1006320219594947, -0.30152637603744714, -1.9609799243212849], [1.6317864231972186, -1.2147664367977435, -1.1287069288319613], [3.280492551617343, -0.695376871636861, -2.431841855490721], [2.3884877116331333, 1.2769614060809018, -0.987748834219553], [1.1325100895164046, 1.9151819587079808, -0.793775261905219], [3.285572867069287, 0.9562100996238483, 0.059969446116429945], [3.196424222201781, 2.073569024157934, -2.080084646834336], [1.23959077527223, -1.0741885454775768, -4.013081746298283], [1.3717373297746525, 1.074587368702989, -3.720168632031241], [-0.7843720924895152, -1.3445493546988831, -2.475221044816154], [-1.1581983763610515, 0.4930222734242889, -3.561495300206362], [-0.2560445759534315, 1.8735740976390405, -1.3916701657561978], [-2.007602477187449, 0.6427130616946731, -1.0543092166280674], [1.5770424436171637, 0.0, 0.0], [2.292718146893912, 1.3915527243580585, 0.0], [3.7823355744197125, 1.318614735245463, -0.48079377647244265], [4.419592497958772, 0.34925073248439586, 0.14976318527135296], [3.8473329977166344, 1.1104432852892008, -1.7795740852228774], [4.368010411051574, 2.468434085335211, -0.20560554420871757], [2.293896732220212, 1.8704189044736097, 1.2405689893126322], [1.649464944000877, 2.235214894314337, -0.8029305726284897], [1.9974224573334838, -0.6906780683055216, 1.0535722235493015], [1.9277183224308936, -0.6529932317206268, -1.1102241252095297], [-0.35014935725347807, 0.7132786644586347, 1.0706086973199291], [-0.42668432219275604, -1.2491488329668503, 0.15362238828850455], [3.8794443758764348, 2.61716360664109, -1.6630105691632127]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0206', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
