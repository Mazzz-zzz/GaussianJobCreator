import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0261'
logfile = 'conf/5009017845242299296281_0261.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863812, -1.3935598728846, 0.08664925740765342], [-0.3466020415139017, -2.4204394252486647, -1.0595513872112599], [-0.7363023803695176, -2.0273092602547957, -2.527598274606108], [0.11232051632789476, -2.8153848177000156, -3.5851242808869475], [-0.4929546906111552, -2.784780389042993, -4.759520822046904], [1.324743353276518, -2.2803342814881122, -3.697604628938583], [0.3402914573864424, -4.612943107426152, -3.0967834407100177], [0.7507765384429437, -5.338265935766086, -4.249102754411298], [1.0389607909617162, -4.618960234585807, -1.865044130567552], [-1.165791972738115, -4.9587360955161985, -2.7923453670381924], [-2.020968275993463, -2.311513884043149, -2.7221562889464304], [-0.5312999584009178, -0.7299509739637076, -2.715388810757052], [0.9809271994195107, -2.5623271837262807, -1.0109564779824434], [-0.9145558447621676, -3.5839764719224823, -0.7664486088494576], [-2.011954187959717, -1.2102165034857002, 0.000966807739598837], [-0.40994706586377416, -1.9496054034114851, 1.2594037813693548], [1.5770424436171655, 0.0, 0.0], [2.2927181468939146, 1.3915527243580565, 0.0], [2.3410798567223186, 2.059852692894932, 1.4165023767064748], [1.1453916410070406, 2.033172555891594, 1.9759511228307947], [3.2022179922846443, 1.4422259961598296, 2.1984030628363236], [2.721555559149498, 3.3152059779287955, 1.2745358845394035], [1.6292323391939714, 2.212255867310577, -0.8090479336198872], [3.54558683009438, 1.2600392214310803, -0.4280914688619854], [1.9974224573334856, -0.6906780683055256, 1.0535722235492975], [1.927718322430895, -0.6529932317206224, -1.1102241252095337], [-0.35014935725347734, 0.5705349971623094, -1.1530217920585795], [-0.42668432219275926, 0.7576153073313049, 1.0049834283127268], [-1.2253125480986944, -5.577404679841883, -2.050948799354474]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0261', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
