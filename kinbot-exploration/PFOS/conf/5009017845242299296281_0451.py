import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0451'
logfile = 'conf/5009017845242299296281_0451.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863848, 0.7718203945763839, 1.1635336229088458], [-2.2709622836291934, 0.7431123812655668, 1.1797556627388963], [-2.9970239643019028, 1.4656877041954859, 2.368157397369676], [-2.2704853996782592, 2.7960258391239927, 2.770559335575194], [-1.2136023415164638, 2.5228411085635236, 3.5154554675534553], [-1.8869755643483324, 3.46097275133602, 1.6844351020282613], [-3.3935951456019624, 3.934431604618535, 3.7525573716297362], [-2.582888557332584, 4.900367948095435, 4.409940681348726], [-4.488149359648491, 4.256070130996564, 2.913595809631594], [-3.899994575451984, 2.8908419465698736, 4.817931893978783], [-4.239914773111564, 1.7601618341371432, 1.9973488528651409], [-3.0230545829338547, 0.6702887991292613, 3.4298036615283585], [-2.6427825570536165, 1.3431159399560042, 0.045483693336761485], [-2.6677893778920354, -0.5232706377056701, 1.1465664052870856], [-0.3710451618282825, 2.061737287215914, 1.0602591291106092], [-0.2455253200204905, 0.3056731502912652, 2.3240905646658523], [1.5770424436171642, 0.0, 0.0], [2.2927181468939186, 1.3915527243580579, 0.0], [2.3410798567223243, 2.0598526928949275, 1.4165023767064726], [1.1453916410070388, 2.0331725558915963, 1.9759511228307927], [3.202217992284642, 1.442225996159828, 2.1984030628363285], [2.721555559149508, 3.31520597792879, 1.2745358845394035], [1.629232339193982, 2.212255867310575, -0.8090479336198879], [3.5455868300943845, 1.260039221431067, -0.42809146886197974], [1.9974224573334796, -0.6906780683055243, 1.0535722235492948], [1.9277183224308931, -0.6529932317206285, -1.1102241252095324], [-0.3501493572534751, -1.2838136616209415, 0.08241309473864945], [-0.42668432219275604, 0.4915335256355461, -1.15860581660123], [-3.8336298501211523, 1.9921374437214996, 4.465954158194528]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0451', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
