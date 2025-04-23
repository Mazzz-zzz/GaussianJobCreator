import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0230'
logfile = 'conf/5009017845242299296281_0230.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863801, -1.3935598728845977, 0.0866492574076508], [-0.34660204151390234, -2.420439425248662, -1.0595513872112663], [1.1624535490467065, -2.799902192493106, -1.2603960465750526], [1.9172819807099781, -2.925959600137114, 0.10851911797417492], [3.0265100569963654, -3.6259410075683354, -0.05316281393413768], [2.2192278579630744, -1.7209079512393501, 0.5834709589555994], [0.8705838155292752, -3.789064502734771, 1.4051054858009535], [-0.09108902253227091, -2.8606418188725313, 1.8905410432487748], [0.5486394924936057, -5.067883464385102, 0.8890331189407575], [1.9831633243341806, -3.9841119286641304, 2.5026782478694978], [1.2315520144708525, -3.9701184155028195, -1.8887358797759306], [1.7685418932081407, -1.8669832707228646, -1.983409881825595], [-1.011797100425895, -3.5349698356523596, -0.7430628651027352], [-0.8215282967123821, -1.946695641113873, -2.2051128650058014], [-2.011954187959715, -1.2102165034856978, 0.0009668077395988351], [-0.4099470658637734, -1.9496054034114874, 1.2594037813693488], [1.5770424436171677, 0.0, 0.0], [2.2927181468939186, 1.3915527243580554, 0.0], [2.3410798567223297, 2.0598526928949283, 1.4165023767064726], [1.145391641007047, 2.0331725558915905, 1.9759511228307989], [3.202217992284658, 1.4422259961598258, 2.1984030628363254], [2.721555559149505, 3.3152059779287937, 1.2745358845394066], [1.6292323391939747, 2.2122558673105743, -0.8090479336198877], [3.545586830094381, 1.260039221431079, -0.42809146886199084], [1.9974224573334864, -0.6906780683055285, 1.0535722235492928], [1.9277183224308942, -0.6529932317206204, -1.110224125209537], [-0.35014935725347474, 0.5705349971623137, -1.1530217920585768], [-0.4266843221927532, 0.7576153073313002, 1.00498342831273], [2.8604463942531058, -3.9884669220981177, 2.0948574098023247]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0230', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
