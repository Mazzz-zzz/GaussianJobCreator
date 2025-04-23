import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0185'
logfile = 'conf/5009017845242299296281_0185.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863831, -1.393559872884596, 0.08664925740765317], [-0.3976197158559577, -2.243200645541632, 1.3824827499919234], [-0.74918335176666, -1.5977521800128298, 2.768620777958028], [-0.9314389667275313, -2.6832583390502456, 3.8857798662301843], [-2.1269824923960123, -3.2369910458976165, 3.7830291451784337], [0.0030362538243395775, -3.623143539359187, 3.7750011996568285], [-0.7791576068935264, -1.9468362525765637, 5.605127483530379], [0.6003697055768613, -1.7340382242896768, 5.877253191826081], [-1.7728499369113317, -0.9436093357118696, 5.712216596842323], [-1.2703030198852923, -3.1945197135590426, 6.431326487582634], [0.24079307197239408, -0.7855357992989355, 3.128124302067116], [-1.8795656402123533, -0.9099255486996073, 2.6708444581582906], [-1.1305489784219354, -3.352053415327125, 1.247992163869577], [0.8878652548597774, -2.574602901845816, 1.382089214148239], [-0.25604457595342883, -2.142008766053287, -0.9267276815498177], [-2.007602477187447, -1.2344150958913065, -0.029451230457632926], [1.577042443617165, 0.0, 0.0], [2.292718146893923, 1.3915527243580512, 0.0], [2.3410798567223328, 2.0598526928949257, 1.4165023767064722], [1.145391641007055, 2.0331725558916047, 1.9759511228307844], [3.202217992284657, 1.4422259961598185, 2.1984030628363236], [2.7215555591495315, 3.31520597792879, 1.2745358845393941], [1.6292323391939947, 2.212255867310572, -0.8090479336198932], [3.54558683009439, 1.2600392214310583, -0.42809146886198507], [1.9974224573334824, -0.6906780683055254, 1.0535722235493004], [1.9277183224308985, -0.6529932317206296, -1.1102241252095302], [-0.3501493572534703, 0.5705349971623088, -1.1530217920585821], [-0.4266843221927538, 0.757615307331308, 1.0049834283127252], [-0.5137090202813567, -3.714793314085684, 6.736028451892618]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0185', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
