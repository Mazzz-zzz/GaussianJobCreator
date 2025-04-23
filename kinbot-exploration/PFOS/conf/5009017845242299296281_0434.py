import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0434'
logfile = 'conf/5009017845242299296281_0434.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863868, 0.7718203945763855, 1.163533622908844], [-0.3466020415139019, 0.2926212946843474, 2.62593772419238], [-0.9873929842445243, 1.0841115196807838, 3.8194659122849317], [-2.4621724959026214, 1.5155496356173197, 3.5056638996836145], [-2.4615979848336345, 2.5854211014867863, 2.729856797671557], [-3.1141435203244794, 0.5282356548563647, 2.89817321888162], [-3.414181650933437, 1.9326408849350447, 5.0681544723256025], [-2.567225413766636, 2.701447264295199, 5.913159141655534], [-4.710125613672749, 2.3364062747563463, 4.664595892847989], [-3.5305767338543412, 0.47731855618511304, 5.658943549681231], [-0.9969507321948189, 0.3041218762439164, 4.896738560267489], [-0.27888546472281034, 2.1786256533403434, 4.065404334434961], [-0.7784534903451367, -0.9699349115591068, 2.691387238785718], [0.973158224737945, 0.3037494144073532, 2.7687950826457906], [-2.011954187959716, 0.6059455318059168, 1.0475948322279929], [-0.40994706586377694, 2.0654783699937846, 1.058705916025071], [1.5770424436171655, 0.0, 0.0], [2.2927181468939195, 1.391552724358057, 0.0], [3.7823355744197125, 1.3186147352454658, -0.48079377647244953], [4.419592497958774, 0.3492507324843962, 0.14976318527135063], [3.847332997716637, 1.1104432852891861, -1.7795740852228794], [4.368010411051578, 2.468434085335211, -0.2056055442087243], [2.2938967322202157, 1.8704189044736084, 1.2405689893126284], [1.649464944000881, 2.235214894314332, -0.8029305726284923], [1.9974224573334851, -0.6906780683055231, 1.0535722235493052], [1.9277183224308976, -0.6529932317206266, -1.1102241252095304], [-0.3501493572534725, -1.2838136616209497, 0.08241309473865455], [-0.4266843221927561, 0.4915335256355397, -1.158605816601228], [-3.463467982016996, -0.1789924643705813, 4.951331617602652]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0434', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
