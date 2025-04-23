import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0216'
logfile = 'conf/5009017845242299296281_0216.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863873, 0.7718203945763832, 1.163533622908847], [-0.3466020415139064, 0.2926212946843425, 2.625937724192379], [-0.7363023803695208, -1.1753096862432062, 3.019500458011146], [0.11232051632789275, -1.6971162941225006, 4.230756914000711], [1.3083470166553104, -2.0726220403760207, 3.812405401367929], [0.2433393664383972, -0.7440121314069538, 5.149085045162456], [-0.6978983773361233, -3.1685933102342845, 5.0673954658972935], [-1.765667469854746, -2.6887668504926183, 5.874901962096378], [-0.8391756380146922, -4.175560712350379, 4.081792017070774], [0.4974230088710358, -3.583678604212362, 6.005187346971094], [-2.0209682759934653, -1.2016995572776057, 3.362907889255021], [-0.5312999584009176, -1.9866202042857686, 1.989850492348295], [0.980927199419505, 0.40564959980990695, 2.7245186729056163], [-0.9145558447621747, 1.1282242700023721, 3.4870389756753206], [-2.0119541879597196, 0.6059455318059146, 1.047594832227991], [-0.4099470658637822, 2.065478369993786, 1.058705916025074], [1.5770424436171635, 0.0, 0.0], [2.292718146893912, 1.3915527243580565, 0.0], [2.3410798567223186, 2.0598526928949306, 1.4165023767064706], [1.1453916410070293, 2.033172555891596, 1.9759511228307922], [3.2022179922846403, 1.4422259961598338, 2.1984030628363262], [2.721555559149497, 3.3152059779287972, 1.2745358845394057], [1.6292323391939716, 2.212255867310582, -0.809047933619891], [3.545586830094379, 1.2600392214310752, -0.4280914688619786], [1.9974224573334816, -0.6906780683055248, 1.0535722235492986], [1.927718322430895, -0.6529932317206224, -1.1102241252095302], [-0.35014935725347657, -1.2838136616209468, 0.08241309473864979], [-0.42668432219276103, 0.4915335256355427, -1.1586058166012294], [0.5355286641010409, -4.546200127573951, 6.094979678498791]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0216', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
