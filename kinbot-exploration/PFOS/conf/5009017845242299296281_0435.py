import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0435'
logfile = 'conf/5009017845242299296281_0435.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863808, 0.7718203945763836, 1.1635336229088498], [-0.3466020415138996, 0.2926212946843405, 2.6259377241923816], [-0.7363023803695162, -1.1753096862432064, 3.0195004580111466], [-2.245266873553545, -1.285267541889903, 3.432441842826776], [-2.9870600258534608, -0.549786506182953, 2.622928094707954], [-2.6496830196655057, -2.551038990399627, 3.376136767970709], [-2.5333983411176324, -0.6864026829068566, 5.187435993583288], [-2.09647104819226, -1.702400504054251, 6.081457232665809], [-2.0946540629194073, 0.6589946099637823, 5.240032137288761], [-4.108201473709184, -0.6739314081950011, 5.1678862954808835], [-0.5266683954723083, -1.970691804482322, 1.9743151008291855], [0.004125934750531798, -1.5779705333152358, 4.044224788693658], [0.9809271994195118, 0.4056495998099024, 2.7245186729056186], [-0.9145558447621631, 1.12822427000237, 3.4870389756753237], [-2.0119541879597147, 0.6059455318059209, 1.0475948322279978], [-0.4099470658637725, 2.065478369993786, 1.0587059160250798], [1.5770424436171668, 0.0, 0.0], [2.2927181468939186, 1.3915527243580543, 0.0], [1.6005215470082543, 2.4407219045638984, -0.9357086002340234], [1.3760692761371698, 1.9145350254105789, -2.1257143081021423], [0.46117091336276017, 2.851946778983255, -0.4188289776134413], [2.404864907693505, 3.4780800111829118, -1.0689303403306805], [3.5400592233304695, 1.2295174124846031, -0.4315210556927465], [2.3088468039522585, 1.896094738758372, 1.2310220414904718], [1.9974224573334824, -0.6906780683055274, 1.0535722235492948], [1.9277183224308942, -0.6529932317206217, -1.110224125209533], [-0.3501493572534787, -1.2838136616209412, 0.08241309473864827], [-0.42668432219275704, 0.4915335256355517, -1.158605816601226], [-4.432311136493154, -0.5821374538838229, 4.260974401668712]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0435', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
