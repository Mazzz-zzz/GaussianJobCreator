import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0464'
logfile = 'conf/5009017845242299296281_0464.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863813, -1.3935598728846004, 0.08664925740765345], [-0.34660204151390006, -2.4204394252486665, -1.0595513872112594], [-0.7363023803695165, -2.0273092602547966, -2.527598274606109], [-2.245266873553545, -2.32994806195571, -2.8292952633496244], [-2.6409584449074375, -1.63499039435252, -3.8814188129748675], [-2.4192802354782867, -3.627374974762621, -3.064460729243292], [-3.3507490824222406, -1.8731712021050504, -1.3832913500441468], [-3.227636012618887, -2.885550614633041, -0.3921495312598397], [-3.154695608522511, -0.49269869212109446, -1.1358737754895982], [-4.743136573073806, -2.0358620541085712, -2.1011572096310283], [-0.5266683954723087, -0.7244611301521506, -2.693826716126078], [0.0041259347505317985, -2.7134161389658433, -3.3886749626211023], [0.9809271994195121, -2.56232718372628, -1.0109564779824443], [-0.9145558447621641, -3.5839764719224823, -0.7664486088494576], [-2.011954187959717, -1.2102165034857038, 0.0009668077396009896], [-0.40994706586377344, -1.9496054034114856, 1.2594037813693526], [1.5770424436171668, 0.0, 0.0], [2.2927181468939093, 1.3915527243580619, 0.0], [3.782335574419716, 1.318614735245461, -0.4807937764724422], [4.4195924979587735, 0.34925073248440297, 0.14976318527134486], [3.8473329977166366, 1.1104432852892059, -1.7795740852228812], [4.368010411051574, 2.468434085335213, -0.20560554420871807], [2.293896732220214, 1.870418904473608, 1.2405689893126348], [1.6494649440008797, 2.2352148943143395, -0.8029305726284804], [1.9974224573334844, -0.6906780683055256, 1.0535722235492992], [1.9277183224308951, -0.6529932317206225, -1.1102241252095313], [-0.3501493572534785, 0.5705349971623092, -1.153021792058579], [-0.42668432219275704, 0.7576153073313022, 1.0049834283127286], [-4.6430901167298995, -1.9268188494403298, -3.0572234587919906]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0464', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
