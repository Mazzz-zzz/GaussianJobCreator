import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0128'
logfile = 'conf/5009017845242299296281_0128.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863838, 0.6217394783082144, -1.250182880316499], [-2.2709622836291943, 0.6501421835576486, -1.2334320314121743], [-2.9699917885127127, 1.3987418462127572, -0.044930604952379785], [-3.056739804078712, 0.5000645662030765, 1.2374262436633887], [-1.9299416078656422, -0.1731733017708804, 1.390979300337152], [-3.2751648795012236, 1.2495736306950789, 2.314193019350682], [-4.462079523721017, -0.7386873768807473, 1.1270757700950764], [-5.680757787577887, -0.05298502843358739, 1.3860479430149752], [-4.2250883308147795, -1.5278386064210359, -0.024632649188188378], [-4.09590600212127, -1.6109083870378045, 2.3863176417499345], [-2.265942241823491, 2.4877087631800117, 0.2508004249777456], [-4.202897212097874, 1.7436062796993288, -0.3934166742568523], [-2.62120547179293, 1.2630031390416223, -2.367697716836083], [-2.71337798074113, -0.6010536646032283, -1.266242407574956], [-0.37104516182828434, -0.11265730320380395, -2.315646431213898], [-0.24552532002049307, 1.8598848945507183, -1.4267659957399774], [1.577042443617164, 0.0, 0.0], [2.2927181468939137, 1.3915527243580594, 0.0], [3.7823355744197125, 1.3186147352454596, -0.4807937764724533], [4.419592497958774, 0.34925073248439586, 0.1497631852713483], [3.847332997716633, 1.110443285289203, -1.7795740852228823], [4.368010411051577, 2.468434085335213, -0.20560554420872307], [2.2938967322202157, 1.870418904473608, 1.2405689893126286], [1.6494649440008766, 2.235214894314339, -0.8029305726284889], [1.9974224573334864, -0.6906780683055218, 1.0535722235492992], [1.9277183224308958, -0.6529932317206212, -1.110224125209533], [-0.3501493572534751, 0.7132786644586355, 1.0706086973199345], [-0.42668432219275365, -1.24914883296685, 0.15362238828850228], [-3.15230135371049, -1.534053737112016, 2.5854809339961435]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0128', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
