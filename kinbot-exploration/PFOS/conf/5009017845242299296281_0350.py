import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0350'
logfile = 'conf/5009017845242299296281_0350.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863821, 0.621739478308218, -1.250182880316502], [-0.3976197158559556, -0.07566485901595241, -2.633910119820669], [-0.7491833517666582, -1.5988198371506641, -2.7680043658220996], [0.3751315092500395, -2.510239815805813, -2.1640052972737625], [-0.10815375808507582, -3.715131603382629, -1.9166733873151727], [1.3965384052652352, -2.6078030699462404, -3.010237774059631], [1.0485429431935351, -1.8142296954322226, -0.5566004683377855], [1.7447219915669392, -2.8529480084722634, 0.12076945022307528], [1.6092624577838988, -0.5495372756666591, -0.8590958914321986], [-0.3233323711935339, -1.5647107042076382, 0.1756843291301464], [-1.8820373528306364, -1.8378201814607653, -2.113395933093814], [-0.8963077255402097, -1.9155699118782799, -4.0480253071706365], [-1.1305489784219276, 0.59523379002861, -3.52695949445047], [0.8878652548597812, 0.09037708117406702, -2.9207161247297284], [-0.25604457595342656, 1.8735740976390463, -1.391670165756194], [-2.007602477187446, 0.6427130616946769, -1.054309216628063], [1.5770424436171644, 0.0, 0.0], [2.292718146893915, 1.3915527243580568, 0.0], [2.3410798567223194, 2.059852692894929, 1.4165023767064742], [1.145391641007034, 2.0331725558915914, 1.9759511228307942], [3.202217992284644, 1.4422259961598294, 2.1984030628363245], [2.72155555914949, 3.315205977928798, 1.2745358845394017], [1.6292323391939716, 2.2122558673105797, -0.8090479336198859], [3.5455868300943756, 1.2600392214310783, -0.428091468861984], [1.9974224573334807, -0.6906780683055282, 1.0535722235492988], [1.9277183224308951, -0.6529932317206247, -1.1102241252095335], [-0.3501493572534758, 0.7132786644586331, 1.070608697319933], [-0.42668432219276087, -1.24914883296685, 0.15362238828849775], [-0.2744489890056814, -0.7646639244666817, 0.7174388416798836]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0350', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
