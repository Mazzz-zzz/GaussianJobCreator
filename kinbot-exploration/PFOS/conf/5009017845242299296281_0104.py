import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0104'
logfile = 'conf/5009017845242299296281_0104.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586382, 0.6217394783082146, -1.2501828803165036], [-2.270962283629193, 0.6501421835576519, -1.2334320314121765], [-3.020318448930583, -0.7271986321270452, -1.1788705129599844], [-2.2823447039971057, -1.7540796922695505, -0.2512417083665389], [-3.1133589271330373, -2.7197724001234005, 0.09997745116306522], [-1.2352156515326025, -2.2783257115331774, -0.8817866308931175], [-1.6404353738143138, -0.9432111831185305, 1.314679108215983], [-1.3547325602404028, -1.9662287051065463, 2.2603217784839083], [-0.7121068992865005, 0.04808607499585326, 0.9135647092176591], [-2.964071167648961, -0.2087858286764197, 1.7495932535841354], [-4.245418844080689, -0.5314114745498671, -0.6995084752088678], [-3.0940580938002604, -1.2477925357801372, -2.397055332165478], [-2.5961406643712635, 1.338449996801186, -0.13556188350402282], [-2.689777026166649, 1.3195442638618275, -2.3005750802147835], [-0.37104516182828245, -0.11265730320380297, -2.3156464312138993], [-0.24552532002048913, 1.8598848945507183, -1.4267659957399803], [1.5770424436171668, 0.0, 0.0], [2.2927181468939195, 1.3915527243580559, 0.0], [3.7823355744197165, 1.3186147352454574, -0.4807937764724494], [4.419592497958777, 0.34925073248439564, 0.14976318527134475], [3.8473329977166366, 1.1104432852891966, -1.779574085222879], [4.36801041105158, 2.468434085335203, -0.20560554420871602], [2.2938967322202193, 1.8704189044736057, 1.2405689893126317], [1.6494649440008824, 2.2352148943143386, -0.802930572628488], [1.9974224573334847, -0.6906780683055261, 1.053572223549299], [1.9277183224308951, -0.6529932317206258, -1.1102241252095346], [-0.3501493572534736, 0.7132786644586335, 1.0706086973199336], [-0.42668432219275704, -1.249148832966851, 0.1536223882885024], [-3.451078465921009, -0.7453971438159701, 2.3905570392927844]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0104', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
