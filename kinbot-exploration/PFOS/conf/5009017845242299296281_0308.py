import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0308'
logfile = 'conf/5009017845242299296281_0308.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863848, 0.7718203945763839, 1.1635336229088458], [-2.2709622836291934, 0.7431123812655668, 1.1797556627388963], [-2.9970239643019028, 1.4656877041954859, 2.368157397369676], [-2.2704853996782592, 2.7960258391239927, 2.770559335575194], [-1.2136023415164638, 2.5228411085635236, 3.5154554675534553], [-1.8869755643483324, 3.46097275133602, 1.6844351020282613], [-3.3935951456019624, 3.934431604618535, 3.7525573716297362], [-2.5828885573325806, 4.900367948095435, 4.409940681348726], [-4.488149359648491, 4.256070130996564, 2.913595809631594], [-3.899994575451984, 2.8908419465698736, 4.817931893978783], [-4.239914773111564, 1.7601618341371432, 1.9973488528651409], [-3.0230545829338547, 0.6702887991292613, 3.4298036615283585], [-2.6427825570536165, 1.3431159399560042, 0.045483693336761485], [-2.6677893778920354, -0.5232706377056701, 1.1465664052870856], [-0.3710451618282825, 2.061737287215914, 1.0602591291106092], [-0.2455253200204905, 0.3056731502912652, 2.3240905646658523], [1.5770424436171642, 0.0, 0.0], [2.2927181468939186, 1.3915527243580579, 0.0], [1.6005215470082543, 2.4407219045638975, -0.9357086002340286], [1.3760692761371678, 1.9145350254105773, -2.125714308102144], [0.46117091336275884, 2.8519467789832533, -0.4188289776134546], [2.4048649076934994, 3.4780800111829073, -1.0689303403306814], [3.5400592233304655, 1.229517412484605, -0.4315210556927418], [2.308846803952248, 1.896094738758372, 1.2310220414904727], [1.99742245733348, -0.6906780683055244, 1.0535722235492946], [1.9277183224308931, -0.6529932317206287, -1.110224125209532], [-0.3501493572534751, -1.2838136616209415, 0.08241309473864945], [-0.42668432219275604, 0.4915335256355461, -1.15860581660123], [-3.366805079684271, 2.9470501231699457, 5.623234314925264]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0308', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
