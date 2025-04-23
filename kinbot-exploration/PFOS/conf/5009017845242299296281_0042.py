import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0042'
logfile = 'conf/5009017845242299296281_0042.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863822, -1.3935598728845964, 0.08664925740764573], [-0.3976197158559575, -2.2432006455416387, 1.382482749991911], [1.0879414097563231, -2.6761952078648323, 1.6417185410974562], [1.4730420899455332, -3.954161328190836, 0.818259240088166], [1.0004737759019138, -5.033376622340047, 1.416987459327409], [0.9773040172258539, -3.878060131557083, -0.4135493536795969], [3.332923574810669, -4.146760933982904, 0.6588553969111794], [3.9198025913470134, -3.858297469635425, 1.921706479777238], [3.561885499597077, -5.340131249902282, -0.06833130382569709], [3.604488890927199, -2.9370329705429943, -0.31240829713028007], [1.8947536492841386, -1.6828188678150067, 1.2794354892121604], [1.2598103705638535, -2.945657386689855, 2.9294611380409665], [-0.7843720924895164, -1.4713296274432173, 2.402024420219272], [-1.158198376361056, -3.3308565421497396, 1.353777836686174], [-0.2560445759534298, -2.142008766053281, -0.9267276815498285], [-2.0076024771874463, -1.2344150958913058, -0.0294512304576384], [1.5770424436171657, 0.0, 0.0], [2.2927181468939186, 1.3915527243580539, 0.0], [1.6005215470082614, 2.4407219045638993, -0.9357086002340307], [1.376069276137172, 1.9145350254105789, -2.1257143081021423], [0.4611709133627675, 2.851946778983258, -0.41882897761344895], [2.40486490769351, 3.4780800111829064, -1.068930340330688], [3.5400592233304664, 1.2295174124845996, -0.43152105569275123], [2.308846803952265, 1.8960947387583702, 1.2310220414904653], [1.9974224573334833, -0.6906780683055292, 1.0535722235492964], [1.9277183224308954, -0.652993231720624, -1.1102241252095366], [-0.35014935725347246, 0.5705349971623128, -1.1530217920585801], [-0.42668432219275704, 0.757615307331302, 1.0049834283127266], [2.8069792691932283, -2.726357062020868, -0.8179342128913991]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0042', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
