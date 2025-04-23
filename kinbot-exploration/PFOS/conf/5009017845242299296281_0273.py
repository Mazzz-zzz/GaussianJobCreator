import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0273'
logfile = 'conf/5009017845242299296281_0273.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863791, -1.3935598728846006, 0.08664925740765346], [-0.3466020415139001, -2.420439425248666, -1.0595513872112616], [-0.9873929842445189, -3.849810268767847, -0.9708648395635693], [-2.462172495902615, -3.7937688120646964, -0.4403274647009681], [-3.0952492719626594, -2.7675473079450414, -0.9812388836715881], [-3.105133356403326, -4.918796415434234, -0.7400177449467298], [-2.5252703668709358, -3.5965120887216355, 1.4248191179809997], [-1.5506217707940506, -2.6325591650291313, 1.8033503358217742], [-3.893028646093177, -3.5252903392798975, 1.784710106578221], [-2.0014308589805956, -5.0275150875359715, 1.8227703928458878], [-0.9969507321948156, -4.392760927004434, -2.1849920094599393], [-0.27888546472280223, -4.610056256946217, -0.14595700608829107], [-0.7784534903451409, -1.8458422644301253, -2.1856818928204627], [0.9731582247379484, -2.5497215866483587, -1.1213428320614909], [-2.0119541879597147, -1.2102165034857026, 0.0009668077396009887], [-0.40994706586376983, -1.949605403411485, 1.2594037813693497], [1.5770424436171684, 0.0, 0.0], [2.2927181468939195, 1.3915527243580577, 0.0], [1.60052154700826, 2.4407219045639024, -0.9357086002340177], [1.3760692761371671, 1.9145350254105789, -2.1257143081021406], [0.46117091336276084, 2.8519467789832573, -0.4188289776134365], [2.4048649076935047, 3.478080011182908, -1.0689303403306796], [3.5400592233304673, 1.2295174124846062, -0.4315210556927468], [2.3088468039522567, 1.896094738758374, 1.2310220414904767], [1.9974224573334876, -0.6906780683055264, 1.0535722235492964], [1.9277183224308947, -0.6529932317206225, -1.110224125209535], [-0.3501493572534782, 0.5705349971623102, -1.153021792058581], [-0.42668432219275415, 0.7576153073313007, 1.0049834283127286], [-1.0488000260677972, -5.002275765178281, 1.9895632568622759]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0273', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
