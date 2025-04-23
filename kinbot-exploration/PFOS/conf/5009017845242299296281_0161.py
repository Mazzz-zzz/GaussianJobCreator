import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0161'
logfile = 'conf/5009017845242299296281_0161.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863829, 0.7718203945763835, 1.1635336229088495], [-0.397619715855956, 2.318865504557592, 1.2514273698287466], [1.0879414097563267, 2.759867566386747, 1.4967937649483787], [1.3412747433709995, 4.235018190162302, 1.0284583955214386], [1.4942413294616161, 4.269507329980803, -0.28375868478629684], [0.3200544491791535, 5.012997571112636, 1.3759598114895844], [2.884934235898717, 4.959924572236463, 1.811586379616129], [2.5770387278690783, 5.317755003815254, 3.153073892111605], [3.9655135290297983, 4.122630042167245, 1.4419337956666904], [2.9826164228846865, 6.273542994671178, 0.9482165513300228], [1.3507067731889262, 2.680176727710471, 2.7981758327967134], [1.9101736356829004, 1.964674530320517, 0.8244605116218737], [-0.7843720924895166, 2.8158789821420864, 0.07319662459687319], [-1.1581983763610533, 2.837834268725446, 2.207717463520174], [-0.2560445759534264, 0.2684346684142451, 2.318397847306015], [-2.007602477187446, 0.5917020341966355, 1.0837604470856936], [1.577042443617165, 0.0, 0.0], [2.292718146893915, 1.3915527243580557, 0.0], [1.6005215470082486, 2.4407219045638997, -0.9357086002340277], [1.3760692761371625, 1.9145350254105762, -2.125714308102146], [0.4611709133627533, 2.8519467789832533, -0.4188289776134513], [2.4048649076934936, 3.4780800111829135, -1.0689303403306836], [3.5400592233304624, 1.2295174124846078, -0.4315210556927457], [2.3088468039522483, 1.8960947387583738, 1.231022041490472], [1.9974224573334867, -0.6906780683055247, 1.0535722235492968], [1.9277183224308965, -0.6529932317206256, -1.110224125209532], [-0.3501493572534734, -1.283813661620947, 0.08241309473865323], [-0.42668432219275737, 0.49153352563554187, -1.158605816601229], [2.535942280220702, 6.151396926678106, 0.09879031259870885]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0161', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
