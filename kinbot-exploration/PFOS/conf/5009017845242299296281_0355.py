import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0355'
logfile = 'conf/5009017845242299296281_0355.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863811, 0.771820394576384, 1.1635336229088502], [-0.3466020415139017, 0.2926212946843379, 2.6259377241923794], [-0.7363023803695176, -1.1753096862432086, 3.0195004580111444], [0.11232051632789476, -1.697116294122508, 4.230756914000707], [-0.4929546906111552, -2.7294757472121183, 4.791450971895395], [1.324743353276518, -2.0620524010676897, 3.823629731358532], [0.3402914573864424, -0.37542157596077735, 5.543317637598385], [0.7507765384429436, -1.0106979607275752, 6.747625289736186], [1.0389607909617167, 0.6943045210423402, 4.932658967505216], [-1.165791972738115, 0.06112602376324484, 5.690564112898983], [-2.020968275993463, -1.2016995572776086, 3.3629078892550193], [-0.5312999584009178, -1.9866202042857697, 1.9898504923482894], [0.9809271994195107, 0.4056495998098998, 2.7245186729056177], [-0.9145558447621676, 1.1282242700023668, 3.487038975675325], [-2.0119541879597165, 0.6059455318059174, 1.047594832227996], [-0.40994706586377416, 2.0654783699937864, 1.058705916025079], [1.5770424436171655, 0.0, 0.0], [2.2927181468939173, 1.3915527243580539, 0.0], [1.6005215470082554, 2.4407219045638984, -0.9357086002340274], [1.376069276137169, 1.9145350254105775, -2.1257143081021432], [0.46117091336276084, 2.8519467789832555, -0.4188289776134456], [2.404864907693506, 3.4780800111829104, -1.0689303403306836], [3.540059223330469, 1.2295174124846038, -0.4315210556927466], [2.3088468039522563, 1.8960947387583738, 1.2310220414904682], [1.9974224573334836, -0.690678068305529, 1.0535722235492968], [1.9277183224308942, -0.6529932317206268, -1.1102241252095328], [-0.3501493572534758, -1.283813661620943, 0.08241309473864725], [-0.42668432219275926, 0.4915335256355472, -1.1586058166012279], [-1.5806666567819585, -0.40674800905356007, 6.4287614257166705]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0355', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
