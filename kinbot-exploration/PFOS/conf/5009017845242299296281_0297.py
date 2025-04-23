import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0297'
logfile = 'conf/5009017845242299296281_0297.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863848, 0.7718203945763861, 1.1635336229088447], [-0.34660204151390644, 0.29262129468434955, 2.6259377241923785], [-0.9873929842445293, 1.0841115196807833, 3.8194659122849295], [-0.18477695124075286, 2.390826822806661, 4.147928663816875], [0.18832441293677163, 2.982413471416188, 3.0266247120211034], [-0.933236841251377, 3.226790268643201, 4.861982548914764], [1.3577185656456119, 2.0276490900586026, 5.153151448047536], [0.96592075388888, 1.7932179546298317, 6.499987564750014], [2.1453375683054725, 1.1327526197350481, 4.3888390867059535], [2.0283879731241963, 3.4498066464179233, 5.062525510417034], [-2.2305068985543404, 1.4246952716906227, 3.4913469454112604], [-1.0036314587762363, 0.32352185066805206, 4.906512125113518], [-0.778453490345142, -0.9699349115591033, 2.691387238785716], [0.9731582247379394, 0.3037494144073555, 2.768795082645788], [-2.011954187959719, 0.6059455318059207, 1.0475948322279887], [-0.40994706586377755, 2.065478369993786, 1.0587059160250682], [1.5770424436171642, 0.0, 0.0], [2.292718146893916, 1.3915527243580557, 0.0], [3.7823355744197142, 1.3186147352454602, -0.48079377647244126], [4.419592497958775, 0.34925073248439575, 0.14976318527135857], [3.8473329977166375, 1.1104432852891943, -1.779574085222876], [4.368010411051578, 2.468434085335205, -0.2056055442087193], [2.293896732220214, 1.870418904473607, 1.2405689893126306], [1.6494649440008806, 2.235214894314333, -0.8029305726284939], [1.9974224573334811, -0.6906780683055226, 1.0535722235493015], [1.9277183224308976, -0.6529932317206273, -1.1102241252095284], [-0.35014935725347474, -1.283813661620943, 0.08241309473865183], [-0.4266843221927567, 0.49153352563554403, -1.158605816601231], [1.8085612480262379, 3.9754672821572963, 5.84439474095548]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0297', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
