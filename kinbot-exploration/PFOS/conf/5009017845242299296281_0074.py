import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0074'
logfile = 'conf/5009017845242299296281_0074.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863868, -1.3935598728845955, 0.08664925740765189], [-0.3466020415139019, -2.4204394252486647, -1.0595513872112636], [-0.9873929842445243, -3.849810268767845, -0.9708648395635706], [-2.4621724959026214, -3.7937688120646906, -0.4403274647009743], [-2.4615979848336345, -3.6568358862206014, 0.8741119545321082], [-3.1141435203244794, -2.77400945954738, -0.9916211131504993], [-3.414181650933437, -5.355470965805208, -0.8603611334166335], [-2.567225413766637, -6.4716696654414765, -0.6170576129641763], [-4.71012561367275, -5.207861678973085, -0.30891075892364966], [-3.5305767338543412, -5.13944815069858, -2.416101779486618], [-0.996950732194819, -4.392760927004432, -2.1849920094599433], [-0.27888546472281034, -4.610056256946215, -0.14595700608828896], [-0.7784534903451367, -1.8458422644301269, -2.185681892820464], [0.973158224737945, -2.549721586648362, -1.1213428320614895], [-2.011954187959716, -1.2102165034856973, 0.0009668077395945351], [-0.40994706586377694, -1.9496054034114818, 1.2594037813693502], [1.5770424436171655, 0.0, 0.0], [2.292718146893919, 1.3915527243580548, 0.0], [1.600521547008265, 2.4407219045638993, -0.9357086002340254], [1.3760692761371738, 1.9145350254105815, -2.1257143081021415], [0.4611709133627613, 2.85194677898326, -0.41882897761344606], [2.4048649076935074, 3.4780800111829056, -1.068930340330683], [3.54005922333047, 1.2295174124846007, -0.4315210556927439], [2.308846803952257, 1.896094738758372, 1.2310220414904716], [1.9974224573334833, -0.690678068305524, 1.0535722235492977], [1.9277183224308962, -0.6529932317206251, -1.1102241252095346], [-0.3501493572534743, 0.5705349971623107, -1.1530217920585861], [-0.42668432219275604, 0.7576153073313042, 1.004983428312724], [-4.380845107241311, -5.469648942366323, -2.7385510268044992]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0074', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
