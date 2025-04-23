import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0232'
logfile = 'conf/5009017845242299296281_0232.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863847, 0.7718203945763868, 1.1635336229088435], [-0.3466020415139065, 0.2926212946843473, 2.6259377241923794], [-0.9873929842445287, 1.0841115196807871, 3.81946591228493], [-2.462172495902626, 1.5155496356173248, 3.505663899683609], [-3.096370318030981, 1.7977976497721508, 4.630222090259049], [-2.4735048074690122, 2.5820291572649725, 2.711061452434815], [-3.4143213947282445, 0.15733163429364252, 2.627985332974563], [-3.008666588967336, 0.13914282215264107, 1.2650677443119265], [-3.39196853699071, -0.9754161325925295, 3.4775379317071953], [-4.852553939297091, 0.7928263778249913, 2.7184077223261807], [-0.9969507321948273, 0.3041218762439206, 4.896738560267487], [-0.2788854647228144, 2.178625653340351, 4.065404334434956], [-0.7784534903451451, -0.9699349115591032, 2.6913872387857154], [0.9731582247379404, 0.30374941440735315, 2.7687950826457897], [-2.01195418795972, 0.6059455318059203, 1.047594832227988], [-0.4099470658637769, 2.0654783699937864, 1.0587059160250671], [1.5770424436171655, 0.0, 0.0], [2.292718146893918, 1.3915527243580519, 0.0], [2.3410798567223288, 2.0598526928949266, 1.4165023767064693], [1.1453916410070428, 2.0331725558915994, 1.975951122830788], [3.2022179922846465, 1.4422259961598227, 2.198403062836327], [2.721555559149519, 3.315205977928792, 1.2745358845393988], [1.6292323391939842, 2.2122558673105726, -0.8090479336198906], [3.5455868300943845, 1.2600392214310605, -0.42809146886197946], [1.9974224573334836, -0.690678068305526, 1.0535722235493064], [1.9277183224309007, -0.6529932317206302, -1.110224125209525], [-0.3501493572534748, -1.2838136616209432, 0.08241309473865299], [-0.42668432219275426, 0.4915335256355418, -1.15860581660123], [-4.901595892137694, 1.4087843918410508, 3.4628222563910196]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0232', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
