import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0122'
logfile = 'conf/5009017845242299296281_0122.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863812, 0.6217394783082113, -1.2501828803165052], [-0.3466020415139017, 2.1278181305643207, -1.5663863369811268], [-0.7363023803695176, 3.2026189464980015, -0.4919021834050462], [0.11232051632789475, 4.51250111182252, -0.6456326331137735], [-0.4929546906111552, 5.514256136255109, -0.031930149848508345], [1.324743353276518, 4.3423866825558015, -0.126025102419963], [0.3402914573864424, 4.988364683386924, -2.4465341968883845], [0.7507765384429437, 6.34896389649365, -2.4985225353249074], [1.0389607909617162, 3.924655713543455, -3.0676148369376763], [-1.165791972738115, 4.897610071752946, -2.898218745860807], [-2.020968275993463, 3.513213441320754, -0.6407516003085993], [-0.5312999584009177, 2.716571178249479, 0.7255383184087529], [0.9809271994195106, 2.156677583916374, -1.7135621949231812], [-0.9145558447621676, 2.4557522019201046, -2.7205903668258746], [-2.011954187959717, 0.6042709716797788, -1.0485616399675963], [-0.40994706586377416, -0.11587296658230846, -2.3181096973944317], [1.5770424436171655, 0.0, 0.0], [2.2927181468939137, 1.3915527243580568, 0.0], [1.600521547008244, 2.4407219045638997, -0.9357086002340324], [1.3760692761371631, 1.9145350254105695, -2.125714308102147], [0.46117091336274707, 2.8519467789832484, -0.41882897761345506], [2.404864907693489, 3.47808001118291, -1.0689303403306907], [3.5400592233304633, 1.2295174124846142, -0.4315210556927461], [2.3088468039522443, 1.896094738758378, 1.2310220414904673], [1.9974224573334836, -0.6906780683055244, 1.0535722235493041], [1.9277183224308958, -0.652993231720628, -1.1102241252095268], [-0.3501493572534768, 0.7132786644586354, 1.0706086973199302], [-0.4266843221927588, -1.2491488329668499, 0.15362238828850225], [-1.2253125480986944, 4.564876102023098, -3.8046997402520604]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0122', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
