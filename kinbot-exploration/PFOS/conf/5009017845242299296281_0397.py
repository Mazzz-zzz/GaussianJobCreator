import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0397'
logfile = 'conf/5009017845242299296281_0397.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863851, -1.393559872884597, 0.08664925740765075], [-0.3976197158559567, -2.2432006455416396, 1.3824827499919141], [-1.1233882121466328, -3.6270084497455466, 1.5236046658601572], [-0.3678111583333635, -4.58148330604694, 2.5123747217585426], [0.6688104135818955, -5.129831292172003, 1.903086492896763], [0.05149838568730977, -3.905908885145356, 3.5785260940644044], [-1.475740760567105, -5.964456490011123, 3.130064368786193], [-2.3402363179501737, -5.430924462791306, 4.125263853817031], [-1.9223669601590485, -6.674239625405266, 1.9889729131192446], [-0.3736644517404418, -6.8409265100765015, 3.8355930664722635], [-2.351535285418464, -3.422611284169413, 1.9914380952175201], [-1.1885301103712402, -4.223582124447544, 0.3402677800844896], [0.9186702629217423, -2.470358014913695, 1.3578310147309456], [-0.6863393184328396, -1.5018460017288953, 2.4451837914307166], [-0.2560445759534336, -2.1420087660532827, -0.9267276815498271], [-2.007602477187449, -1.2344150958913087, -0.02945123045763408], [1.5770424436171657, 0.0, 0.0], [2.2927181468939213, 1.3915527243580526, 0.0], [1.6005215470082612, 2.440721904563901, -0.9357086002340143], [1.3760692761371724, 1.9145350254105862, -2.1257143081021352], [0.46117091336276617, 2.85194677898326, -0.41882897761344073], [2.4048649076935167, 3.4780800111829087, -1.06893034033068], [3.540059223330471, 1.2295174124845947, -0.4315210556927443], [2.3088468039522634, 1.8960947387583669, 1.2310220414904738], [1.997422457333485, -0.6906780683055307, 1.0535722235492928], [1.9277183224308925, -0.652993231720627, -1.11022412520954], [-0.35014935725348195, 0.570534997162314, -1.15302179205858], [-0.4266843221927512, 0.7576153073313039, 1.004983428312731], [0.49451253630617636, -6.66695704146053, 3.445765247926971]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0397', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
