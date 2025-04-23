import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0058'
logfile = 'conf/5009017845242299296281_0058.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863827, -1.393559872884598, 0.08664925740765454], [-0.3976197158559555, -2.2432006455416333, 1.382482749991924], [-0.7491833517666558, -1.5977521800128336, 2.7686207779580285], [-2.0819874597153203, -0.7736518441984128, 2.704769346120318], [-1.8542445281745772, 0.4066120608580374, 2.1555770167573542], [-2.9956443258697933, -1.425410821491696, 1.9910877595570866], [-2.8081615356386624, -0.48876666843532235, 4.411571906231566], [-1.7481451191270514, -0.14103702073414648, 5.2935471307393085], [-3.9802865149534536, 0.2861638277002311, 4.235610175524548], [-3.2557013693023076, -1.9673815655675635, 4.718082449397589], [-0.895313736421605, -2.5657760365450364, 3.6689358763380975], [0.2259041467951681, -0.7802008545701373, 3.144445207584873], [-1.1305489784219298, -3.3520534153271275, 1.247992163869578], [0.8878652548597812, -2.574602901845815, 1.3820892141482384], [-0.2560445759534278, -2.1420087660532916, -0.9267276815498173], [-2.0076024771874463, -1.2344150958913105, -0.02945123045762864], [1.5770424436171662, 0.0, 0.0], [2.2927181468939186, 1.391552724358057, 0.0], [2.3410798567223337, 2.059852692894922, 1.416502376706476], [1.1453916410070437, 2.0331725558915945, 1.9759511228307955], [3.202217992284649, 1.4422259961598183, 2.19840306283633], [2.7215555591495058, 3.3152059779287995, 1.2745358845394015], [1.629232339193983, 2.2122558673105788, -0.8090479336198856], [3.545586830094382, 1.260039221431071, -0.42809146886198013], [1.997422457333484, -0.6906780683055243, 1.0535722235492995], [1.9277183224308958, -0.6529932317206257, -1.1102241252095322], [-0.35014935725347457, 0.5705349971623093, -1.1530217920585808], [-0.4266843221927564, 0.7576153073313067, 1.0049834283127252], [-4.080681532117272, -1.9708870569634338, 5.223411875799684]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0058', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
