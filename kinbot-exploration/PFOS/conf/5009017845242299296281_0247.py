import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0247'
logfile = 'conf/5009017845242299296281_0247.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586384, 0.771820394576383, 1.1635336229088467], [-0.34660204151390517, 0.29262129468434245, 2.6259377241923785], [1.1624535490467045, 0.3084161010830979, 3.0549844500983028], [1.322955526797157, 0.3508819546540591, 4.614478254324753], [1.1329590874195947, 1.583343225134952, 5.052019695636847], [0.4456504601938104, -0.46458122876160457, 5.192626727607083], [3.0339402189610625, -0.19904667675480195, 5.154638614451783], [3.9877413956025296, 0.39364038653776945, 4.28204686096509], [3.074431482668305, -0.089914691196727, 6.565952696766402], [2.904675852559841, -1.7276879343224567, 4.798142809357428], [1.7500514572469128, -0.7950580504033853, 2.601040568339858], [1.763328214900581, 1.3769084895275183, 2.5473394429719383], [-1.011797100425901, 1.1239736000383729, 3.4329051118380063], [-0.821528296712382, -0.9363359387499641, 2.7884443111439516], [-2.011954187959717, 0.6059455318059163, 1.0475948322279898], [-0.40994706586377744, 2.0654783699937864, 1.0587059160250731], [1.5770424436171686, 0.0, 0.0], [2.292718146893913, 1.3915527243580559, 0.0], [2.3410798567223146, 2.059852692894934, 1.4165023767064742], [1.1453916410070264, 2.0331725558915887, 1.9759511228307924], [3.2022179922846363, 1.4422259961598323, 2.1984030628363307], [2.721555559149488, 3.315205977928801, 1.2745358845394068], [1.6292323391939711, 2.212255867310578, -0.8090479336198874], [3.545586830094378, 1.2600392214310865, -0.42809146886197935], [1.9974224573334838, -0.6906780683055241, 1.0535722235493012], [1.927718322430898, -0.652993231720627, -1.1102241252095288], [-0.3501493572534734, -1.2838136616209452, 0.08241309473865081], [-0.4266843221927547, 0.4915335256355483, -1.158605816601229], [1.975000503227351, -1.9953930458125129, 4.798466571006524]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0247', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
