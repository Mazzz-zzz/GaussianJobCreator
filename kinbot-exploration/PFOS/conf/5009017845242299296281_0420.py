import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0420'
logfile = 'conf/5009017845242299296281_0420.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863829, 0.6217394783082103, -1.250182880316506], [-0.34660204151390533, 2.1278181305643202, -1.566386336981128], [-0.736302380369522, 3.2026189464980006, -0.4919021834050475], [0.11232051632789168, 4.512501111822519, -0.6456326331137754], [-0.4929546906111587, 5.514256136255106, -0.031930149848510773], [1.3247433532765127, 4.342386682555802, -0.12602510241996304], [0.3402914573864387, 4.988364683386923, -2.4465341968883867], [1.3354555142780866, 4.14156110448163, -3.007785661415776], [-0.954546671528242, 5.161298161115175, -2.9932013751948423], [0.960169919753806, 6.419477350845224, -2.226966770473475], [-2.020968275993467, 3.513213441320754, -0.6407516003086009], [-0.5312999584009207, 2.7165711782494784, 0.7255383184087527], [0.9809271994195091, 2.156677583916374, -1.7135621949231812], [-0.9145558447621691, 2.455752201920102, -2.720590366825877], [-2.011954187959719, 0.6042709716797774, -1.048561639967597], [-0.40994706586377266, -0.11587296658231003, -2.318109697394432], [1.577042443617163, 0.0, 0.0], [2.2927181468939106, 1.391552724358057, 0.0], [3.7823355744197116, 1.3186147352454605, -0.4807937764724469], [4.41959249795877, 0.3492507324843944, 0.14976318527135513], [3.8473329977166344, 1.1104432852891957, -1.779574085222877], [4.368010411051576, 2.468434085335206, -0.20560554420871852], [2.293896732220212, 1.8704189044736086, 1.240568989312631], [1.6494649440008802, 2.235214894314337, -0.8029305726284902], [1.9974224573334831, -0.6906780683055229, 1.0535722235493048], [1.9277183224308947, -0.6529932317206284, -1.1102241252095268], [-0.3501493572534797, 0.7132786644586351, 1.0706086973199298], [-0.4266843221927602, -1.2491488329668512, 0.15362238828850128], [1.9261217123109915, 6.375740401375758, -2.2583769663127535]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0420', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
