import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0333'
logfile = 'conf/5009017845242299296281_0333.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863808, 0.7718203945763853, 1.1635336229088502], [-0.34660204151390184, 0.29262129468434267, 2.6259377241923803], [-0.7363023803695191, -1.175309686243206, 3.0195004580111457], [-0.5090926255515547, -2.175074851556421, 1.832727311810156], [-0.44469526072129173, -3.4119591218724095, 2.2936744067405836], [-1.5014279159495572, -2.081639779175608, 0.9521362607354197], [1.0867274347625948, -1.814917428111401, 0.913320820343304], [2.1002584267399214, -1.5401402239983664, 1.8723220842691006], [1.2029631916150896, -2.790644773434742, -0.10638608267598974], [0.6562002816258952, -0.4619118152849617, 0.23173090559898477], [0.01837074119013481, -1.5602715505038287, 4.044812310652281], [-2.0169579203546344, -1.2240697706113262, 3.3629733733434097], [0.9809271994195093, 0.4056495998099023, 2.7245186729056177], [-0.9145558447621683, 1.1282242700023728, 3.4870389756753224], [-2.011954187959715, 0.6059455318059227, 1.0475948322279944], [-0.4099470658637725, 2.0654783699937873, 1.0587059160250758], [1.577042443617166, 0.0, 0.0], [2.2927181468939177, 1.391552724358053, 0.0], [2.341079856722328, 2.0598526928949257, 1.4165023767064735], [1.1453916410070413, 2.033172555891596, 1.97595112283079], [3.20221799228465, 1.442225996159825, 2.1984030628363245], [2.721555559149509, 3.3152059779287946, 1.274535884539401], [1.629232339193979, 2.2122558673105748, -0.8090479336198871], [3.545586830094382, 1.2600392214310714, -0.42809146886198374], [1.9974224573334827, -0.6906780683055291, 1.0535722235492975], [1.9277183224308956, -0.6529932317206261, -1.1102241252095326], [-0.3501493572534776, -1.2838136616209415, 0.08241309473864945], [-0.42668432219275615, 0.49153352563554903, -1.158605816601227], [1.047190797369724, -0.3906483029489911, -0.6503180245085912]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0333', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
