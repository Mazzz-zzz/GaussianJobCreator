import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0202'
logfile = 'conf/5009017845242299296281_0202.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586385, 0.7718203945763881, 1.1635336229088498], [-0.39761971585595557, 2.3188655045575928, 1.251427369828747], [1.0879414097563278, 2.7598675663867454, 1.4967937649483825], [2.1006320219594947, 1.8490216187922408, 0.719360460601157], [2.2762575910635956, 0.7136703779757698, 1.3725336095673066], [1.6464483143989894, 1.5907192761912547, -0.5037082619759717], [3.7756073465590045, 2.673466043561305, 0.5285401693191063], [4.727954801297502, 1.6728239946283732, 0.19063363909158315], [3.5683827015011587, 3.856867461195935, -0.22111028447185457], [3.9951505101317784, 3.0929209439357512, 2.0306723258517216], [1.23959077527223, 4.0125250124967184, 1.0762663043113099], [1.3717373297746525, 2.6844668573495607, 2.7907042758982867], [-0.7843720924895146, 2.8158789821420935, 0.07319662459687462], [-1.158198376361052, 2.8378342687254494, 2.2077174635201766], [-0.2560445759534302, 0.2684346684142453, 2.318397847306017], [-2.007602477187449, 0.5917020341966392, 1.0837604470856936], [1.577042443617164, 0.0, 0.0], [2.2927181468939146, 1.3915527243580534, 0.0], [2.341079856722322, 2.05985269289493, 1.4165023767064744], [1.1453916410070346, 2.0331725558915914, 1.9759511228307927], [3.2022179922846403, 1.442225996159828, 2.1984030628363245], [2.7215555591494986, 3.315205977928796, 1.2745358845394057], [1.6292323391939718, 2.212255867310576, -0.8090479336198855], [3.5455868300943783, 1.260039221431075, -0.4280914688619814], [1.9974224573334805, -0.6906780683055288, 1.0535722235493021], [1.9277183224308914, -0.6529932317206287, -1.1102241252095304], [-0.35014935725347923, -1.283813661620942, 0.0824130947386506], [-0.42668432219276153, 0.4915335256355458, -1.1586058166012259], [4.52729672489061, 2.4287368261245463, 2.490718247973855]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0202', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
