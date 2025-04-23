import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0222'
logfile = 'conf/5009017845242299296281_0222.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863813, 0.6217394783082112, -1.2501828803165023], [-0.346602041513904, 2.1278181305643216, -1.5663863369811248], [1.1624535490467054, 2.4914860914100143, -1.7945884035232624], [1.9172819807099761, 1.3689994871066287, -2.5882149031527932], [1.1464470566547256, 0.9090692451221105, -3.5581528575817325], [3.0445827205004465, 1.8471459380411672, -3.1073189910576895], [2.381891830767319, -0.0744894972914171, -1.4827102354710175], [1.2743295608764538, -0.3626569882110682, -0.6384968007868369], [3.0228949468799033, -1.0326122514290719, -2.30518975788215], [3.4947646482606083, 0.6418311068730521, -0.6289295177553924], [1.2315520144708516, 3.6207524606765054, -2.4938554639699166], [1.7685418932081398, 2.651174979139487, -0.625149999973776], [-1.0117971004258972, 2.4109962356139842, -2.6898422467352874], [-0.821528296712383, 2.883031579863844, -0.5833314461381615], [-2.0119541879597174, 0.6042709716797763, -1.0485616399675939], [-0.40994706586377505, -0.11587296658230632, -2.318109697394431], [1.5770424436171668, 0.0, 0.0], [2.292718146893914, 1.3915527243580597, 0.0], [2.3410798567223168, 2.059852692894932, 1.4165023767064742], [1.1453916410070277, 2.0331725558915927, 1.975951122830791], [3.2022179922846434, 1.4422259961598336, 2.1984030628363294], [2.7215555591494907, 3.3152059779288, 1.2745358845394077], [1.629232339193972, 2.212255867310579, -0.8090479336198868], [3.545586830094379, 1.2600392214310838, -0.42809146886198035], [1.9974224573334847, -0.690678068305522, 1.0535722235493006], [1.9277183224308956, -0.652993231720624, -1.1102241252095302], [-0.3501493572534756, 0.7132786644586351, 1.0706086973199318], [-0.42668432219275343, -1.249148832966852, 0.15362238828850253], [3.1203561695082196, 0.9771045933662271, 0.19773387625031144]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0222', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
