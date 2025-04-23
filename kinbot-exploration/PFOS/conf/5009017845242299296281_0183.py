import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0183'
logfile = 'conf/5009017845242299296281_0183.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586388, -1.393559872884595, 0.08664925740765311], [-0.3466020415139043, -2.4204394252486656, -1.0595513872112614], [-0.9873929842445259, -3.8498102687678455, -0.9708648395635672], [-2.4621724959026237, -3.7937688120646906, -0.4403274647009709], [-2.461597984833636, -3.6568358862205987, 0.8741119545321144], [-3.114143520324482, -2.774009459547378, -0.9916211131504972], [-3.4141816509334393, -5.355470965805207, -0.8603611334166283], [-4.567822843373353, -5.409857415587863, -0.030550841243862707], [-3.447907588227803, -5.449055223219134, -2.272969118515789], [-2.3659963988850508, -6.4083881622447, -0.3376607559160367], [-0.9969507321948222, -4.392760927004433, -2.184992009459939], [-0.27888546472281345, -4.610056256946216, -0.14595700608828488], [-0.7784534903451381, -1.8458422644301269, -2.185681892820464], [0.9731582247379437, -2.549721586648363, -1.1213428320614873], [-2.011954187959716, -1.210216503485696, 0.000966807739594534], [-0.40994706586377944, -1.949605403411481, 1.2594037813693497], [1.5770424436171655, 0.0, 0.0], [2.2927181468939204, 1.3915527243580526, 0.0], [2.34107985672233, 2.05985269289493, 1.4165023767064648], [1.1453916410070477, 2.0331725558916007, 1.9759511228307893], [3.2022179922846536, 1.4422259961598247, 2.198403062836327], [2.7215555591495137, 3.315205977928791, 1.2745358845393986], [1.6292323391939836, 2.2122558673105734, -0.8090479336198915], [3.5455868300943854, 1.2600392214310665, -0.42809146886198446], [1.9974224573334816, -0.6906780683055243, 1.053572223549299], [1.9277183224308956, -0.6529932317206266, -1.1102241252095337], [-0.35014935725347307, 0.5705349971623098, -1.153021792058587], [-0.42668432219275587, 0.7576153073313051, 1.0049834283127232], [-2.3647343815181263, -7.19168945360504, -0.9054634337552538]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0183', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
