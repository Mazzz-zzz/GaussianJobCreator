import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0453'
logfile = 'conf/5009017845242299296281_0453.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586381, 0.7718203945763842, 1.1635336229088462], [-0.346602041513899, 0.2926212946843428, 2.6259377241923816], [-0.7363023803695148, -1.1753096862432064, 3.0195004580111466], [-0.5090926255515525, -2.175074851556422, 1.8327273118101537], [-1.5087734085102666, -2.086021433268666, 0.9729542231691364], [0.6350825859976155, -1.9037535085902861, 1.2113193502579416], [-0.4118545896983535, -3.952415117627549, 2.4270922452480397], [-0.6081215815946575, -4.807150477617645, 1.3075665688633635], [0.6988372542698227, -4.033254778941512, 3.301840979770827], [-1.716126564814229, -3.968783177130525, 3.3097866042219977], [0.018370741190139622, -1.5602715505038285, 4.04481231065228], [-2.0169579203546295, -1.2240697706113273, 3.3629733733434124], [0.9809271994195142, 0.40564959980990223, 2.7245186729056172], [-0.9145558447621636, 1.1282242700023704, 3.4870389756753255], [-2.0119541879597134, 0.6059455318059226, 1.0475948322279962], [-0.4099470658637751, 2.0654783699937855, 1.0587059160250716], [1.5770424436171673, 0.0, 0.0], [2.2927181468939173, 1.3915527243580552, 0.0], [1.600521547008248, 2.4407219045638993, -0.935708600234024], [1.3760692761371613, 1.9145350254105749, -2.1257143081021437], [0.4611709133627526, 2.851946778983251, -0.4188289776134464], [2.4048649076934936, 3.4780800111829118, -1.0689303403306867], [3.5400592233304655, 1.2295174124846076, -0.4315210556927453], [2.3088468039522563, 1.8960947387583746, 1.2310220414904696], [1.997422457333485, -0.6906780683055282, 1.0535722235492984], [1.9277183224308962, -0.6529932317206252, -1.1102241252095333], [-0.35014935725347796, -1.2838136616209415, 0.08241309473864945], [-0.4266843221927597, 0.4915335256355513, -1.1586058166012267], [-1.5850071069501621, -4.521503251813428, 4.0929025557310705]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0453', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
