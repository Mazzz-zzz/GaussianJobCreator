import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0134'
logfile = 'conf/5009017845242299296281_0134.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863843, -1.3935598728845984, 0.08664925740765084], [-0.3976197158559566, -2.2432006455416373, 1.3824827499919157], [1.0879414097563251, -2.676195207864829, 1.6417185410974604], [2.100632021959493, -1.5474952427548003, 1.2416194637201279], [3.2636602341713647, -1.7574463535228209, 1.833048116867345], [2.2719309703990027, -1.525272508610281, -0.07711586260721295], [1.4886616557511099, 0.15148205615249102, 1.7521405110780743], [2.592881324393371, 1.0473655484210296, 1.7292778893058112], [0.2718090320333246, 0.37441003423396235, 1.062989163123471], [1.134464322195662, -0.17538060765754604, 3.251556289879345], [1.2395907752722295, -2.9383364670191545, 2.9368154419869685], [1.3717373297746493, -3.7590542260525557, 0.9294643561329435], [-0.7843720924895169, -1.471329627443215, 2.402024420219278], [-1.1581983763610537, -3.3308565421497427, 1.3537778366861783], [-0.25604457595343105, -2.1420087660532845, -0.9267276815498257], [-2.0076024771874494, -1.2344150958913092, -0.02945123045763354], [1.5770424436171655, 0.0, 0.0], [2.292718146893915, 1.3915527243580563, 0.0], [3.7823355744197116, 1.3186147352454596, -0.48079377647244936], [4.419592497958778, 0.3492507324843903, 0.1497631852713478], [3.8473329977166366, 1.1104432852892, -1.7795740852228805], [4.368010411051579, 2.4684340853352063, -0.20560554420871774], [2.2938967322202197, 1.8704189044736004, 1.2405689893126326], [1.6494649440008837, 2.2352148943143373, -0.8029305726284834], [1.9974224573334816, -0.6906780683055272, 1.0535722235492966], [1.9277183224308936, -0.6529932317206237, -1.1102241252095344], [-0.350149357253477, 0.5705349971623103, -1.1530217920585775], [-0.4266843221927596, 0.7576153073313039, 1.004983428312729], [0.9564676515766036, -1.1205116273896842, 3.356470194985506]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0134', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
