import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0158'
logfile = 'conf/5009017845242299296281_0158.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863841, -1.393559872884598, 0.08664925740765081], [-0.3976197158559555, -2.2432006455416396, 1.3824827499919141], [1.0879414097563265, -2.6761952078648297, 1.6417185410974577], [2.1006320219594947, -1.547495242754801, 1.2416194637201239], [1.6317864231972183, -0.3701056553971297, 1.616372058347528], [3.280492551617343, -1.7583483890228284, 1.818134963787024], [2.3884877116331333, -1.493896286033045, -0.612006600208588], [2.973467943634744, -0.23858511692773934, -0.9351520855504238], [2.941815915653076, -2.744838114663385, -0.9785160975849094], [0.8833722872031099, -1.484035903803817, -1.0757661950571806], [1.23959077527223, -2.938336467019159, 2.9368154419869628], [1.3717373297746511, -3.7590542260525557, 0.92946435613294], [-0.7843720924895157, -1.47132962744322, 2.402024420219276], [-1.158198376361052, -3.330856542149744, 1.3537778366861755], [-0.25604457595343105, -2.1420087660532836, -0.9267276815498274], [-2.0076024771874494, -1.2344150958913127, -0.029451230457634172], [1.5770424436171633, 0.0, 0.0], [2.2927181468939106, 1.3915527243580565, 0.0], [3.78233557441971, 1.3186147352454638, -0.4807937764724471], [4.419592497958776, 0.34925073248439387, 0.1497631852713484], [3.847332997716635, 1.110443285289204, -1.779574085222881], [4.368010411051576, 2.4684340853352102, -0.20560554420871385], [2.293896732220217, 1.8704189044736008, 1.240568989312635], [1.6494649440008793, 2.235214894314341, -0.8029305726284803], [1.997422457333482, -0.6906780683055274, 1.053572223549295], [1.9277183224308918, -0.6529932317206204, -1.1102241252095346], [-0.3501493572534775, 0.5705349971623132, -1.1530217920585757], [-0.4266843221927602, 0.757615307331301, 1.0049834283127288], [0.5863374461071714, -0.5769675875324622, -1.2337527751447424]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0158', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
