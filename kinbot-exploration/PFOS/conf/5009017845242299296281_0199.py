import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0199'
logfile = 'conf/5009017845242299296281_0199.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863848, 0.7718203945763861, 1.1635336229088447], [-0.34660204151390644, 0.29262129468434955, 2.6259377241923785], [-0.9873929842445293, 1.0841115196807833, 3.8194659122849295], [-0.18477695124075286, 2.390826822806661, 4.147928663816875], [0.18832441293677163, 2.982413471416188, 3.0266247120211034], [-0.933236841251377, 3.226790268643201, 4.861982548914764], [1.3577185656456119, 2.0276490900586026, 5.153151448047536], [0.9659207538888785, 1.7932179546298317, 6.499987564750014], [2.1453375683054725, 1.1327526197350481, 4.3888390867059535], [2.0283879731241963, 3.4498066464179233, 5.062525510417034], [-2.2305068985543404, 1.4246952716906227, 3.4913469454112604], [-1.0036314587762363, 0.32352185066805206, 4.906512125113518], [-0.778453490345142, -0.9699349115591033, 2.691387238785716], [0.9731582247379394, 0.3037494144073555, 2.768795082645788], [-2.011954187959719, 0.6059455318059207, 1.0475948322279887], [-0.40994706586377755, 2.065478369993786, 1.0587059160250682], [1.5770424436171642, 0.0, 0.0], [2.292718146893916, 1.3915527243580557, 0.0], [2.3410798567223217, 2.05985269289493, 1.416502376706472], [1.1453916410070348, 2.0331725558916003, 1.9759511228307884], [3.202217992284644, 1.4422259961598254, 2.198403062836327], [2.721555559149505, 3.3152059779287946, 1.2745358845394026], [1.6292323391939791, 2.212255867310575, -0.8090479336198886], [3.5455868300943827, 1.2600392214310678, -0.428091468861978], [1.9974224573334818, -0.6906780683055231, 1.0535722235493012], [1.9277183224308971, -0.6529932317206272, -1.1102241252095286], [-0.35014935725347474, -1.283813661620943, 0.08241309473865183], [-0.4266843221927567, 0.49153352563554403, -1.158605816601231], [1.7153424819957388, 3.919535957241338, 4.276828824835692]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0199', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
