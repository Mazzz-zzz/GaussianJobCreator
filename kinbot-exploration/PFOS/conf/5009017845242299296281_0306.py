import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0306'
logfile = 'conf/5009017845242299296281_0306.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863837, 0.621739478308212, -1.2501828803164998], [-0.397619715855957, -0.07566485901595642, -2.6339101198206665], [1.0879414097563251, -0.0836723585219196, -3.1385123060458486], [1.3412747433709995, -1.2268379978241974, -4.18186253593047], [2.433578518280429, -0.9652098327338772, -4.87822817955965], [1.4798103010878627, -2.396013338547869, -3.5632051489667536], [-0.08698699647460374, -1.400399692621987, -5.386689524218467], [-0.4971995188277743, -0.09704445132522034, -5.780925252303704], [0.2629292866236743, -2.430192643116133, -6.293575319444772], [-1.155974242147192, -1.986246685126051, -4.389397089787981], [1.3507067731889257, 1.0832029916024, -3.7201890492274705], [1.9101736356828996, -0.26833351767859714, -2.113688309236765], [-0.7843720924895171, -1.3445493546988785, -2.4752210448161533], [-1.1581983763610553, 0.49302227342429444, -3.561495300206356], [-0.25604457595342917, 1.8735740976390423, -1.3916701657561925], [-2.0076024771874477, 0.6427130616946739, -1.0543092166280603], [1.577042443617163, 0.0, 0.0], [2.2927181468939137, 1.3915527243580565, 0.0], [2.3410798567223186, 2.0598526928949314, 1.416502376706474], [1.1453916410070366, 2.0331725558915967, 1.9759511228307929], [3.2022179922846434, 1.4422259961598298, 2.198403062836325], [2.7215555591494964, 3.3152059779287995, 1.2745358845393975], [1.6292323391939718, 2.212255867310578, -0.8090479336198898], [3.5455868300943774, 1.2600392214310756, -0.4280914688619855], [1.997422457333484, -0.6906780683055287, 1.0535722235492992], [1.9277183224308954, -0.652993231720625, -1.1102241252095322], [-0.35014935725347457, 0.7132786644586319, 1.0706086973199331], [-0.4266843221927553, -1.2491488329668523, 0.15362238828850255], [-0.7178336349735206, -2.4262356317387708, -3.6475042930173704]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0306', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
