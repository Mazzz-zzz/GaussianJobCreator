import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0141'
logfile = 'conf/5009017845242299296281_0141.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863845, 0.7718203945763858, 1.1635336229088464], [-0.3976197158559566, 2.318865504557592, 1.2514273698287466], [1.0879414097563251, 2.759867566386746, 1.4967937649483813], [2.1006320219594943, 1.8490216187922413, 0.7193604606011562], [1.6317864231972183, 1.5848720921948654, -0.48766512951556645], [3.280492551617343, 2.4537252606596796, 0.6137068917036919], [2.388487711633133, 0.21693487995214347, 1.599755434428134], [3.264157839292738, 0.4510836814201246, 2.6955689103927765], [1.1272830289191942, -0.416931536562895, 1.7134782815733656], [3.1903514896851033, -0.5177372628732585, 0.4605390573453791], [1.239590775272229, 4.012525012496718, 1.0762663043113117], [1.3717373297746491, 2.68446685734956, 2.7907042758982863], [-0.7843720924895158, 2.8158789821420913, 0.07319662459687332], [-1.1581983763610537, 2.837834268725447, 2.2077174635201753], [-0.2560445759534325, 0.2684346684142449, 2.3183978473060134], [-2.00760247718745, 0.5917020341966386, 1.0837604470856947], [1.5770424436171646, 0.0, 0.0], [2.292718146893914, 1.3915527243580545, 0.0], [2.3410798567223203, 2.059852692894933, 1.4165023767064746], [1.1453916410070297, 2.0331725558915927, 1.975951122830791], [3.2022179922846363, 1.4422259961598307, 2.198403062836327], [2.7215555591494955, 3.3152059779287972, 1.2745358845394053], [1.6292323391939718, 2.212255867310576, -0.8090479336198874], [3.5455868300943783, 1.2600392214310767, -0.4280914688619804], [1.9974224573334798, -0.6906780683055278, 1.053572223549301], [1.9277183224308945, -0.6529932317206266, -1.1102241252095306], [-0.3501493572534773, -1.2838136616209417, 0.08241309473864603], [-0.42668432219276015, 0.4915335256355453, -1.1586058166012283], [2.973829451338112, -0.13692460364443984, -0.402051408737538]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0141', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
