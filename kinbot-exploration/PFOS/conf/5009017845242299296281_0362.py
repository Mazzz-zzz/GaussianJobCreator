import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0362'
logfile = 'conf/5009017845242299296281_0362.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863828, 0.7718203945763844, 1.1635336229088464], [-0.39761971585595524, 2.3188655045575914, 1.251427369828745], [1.0879414097563267, 2.7598675663867462, 1.4967937649483798], [1.3412747433710015, 4.235018190162302, 1.0284583955214386], [0.3135944378492595, 4.995291633511231, 1.3636589683999778], [2.4483457429192916, 4.71481705577196, 1.5880000749983498], [1.563806489688203, 4.347671144906445, -0.8315018902866815], [0.5864285477808986, 3.5183315255715337, -1.447423283229367], [1.7514335240133216, 5.716749686410805, -1.1409801848901342], [2.9662550528587026, 3.638641198043492, -0.9361784966585324], [1.3507067731889266, 2.6801767277104696, 2.798175832796712], [1.9101736356829009, 1.964674530320516, 0.8244605116218733], [-0.7843720924895167, 2.815878982142089, 0.07319662459687076], [-1.1581983763610535, 2.8378342687254463, 2.20771746352017], [-0.2560445759534288, 0.26843466841424485, 2.318397847306013], [-2.007602477187445, 0.5917020341966361, 1.0837604470856927], [1.5770424436171642, 0.0, 0.0], [2.292718146893915, 1.391552724358054, 0.0], [2.341079856722322, 2.0598526928949292, 1.4165023767064755], [1.145391641007034, 2.033172555891593, 1.9759511228307938], [3.2022179922846394, 1.4422259961598285, 2.1984030628363285], [2.7215555591495035, 3.3152059779287955, 1.2745358845394064], [1.6292323391939754, 2.2122558673105757, -0.8090479336198858], [3.54558683009438, 1.260039221431073, -0.4280914688619788], [1.997422457333484, -0.6906780683055265, 1.0535722235493], [1.9277183224308958, -0.6529932317206282, -1.1102241252095306], [-0.35014935725347723, -1.2838136616209455, 0.08241309473865083], [-0.4266843221927583, 0.4915335256355428, -1.1586058166012296], [2.853846072602825, 2.700499654112327, -1.1440656088234888]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0362', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
