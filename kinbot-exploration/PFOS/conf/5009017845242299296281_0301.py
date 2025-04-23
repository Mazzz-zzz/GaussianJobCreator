import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0301'
logfile = 'conf/5009017845242299296281_0301.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586382, 0.7718203945763844, 1.1635336229088464], [-0.34660204151390284, 0.2926212946843425, 2.625937724192379], [1.1624535490467065, 0.30841610108309514, 3.054984450098302], [1.9243606525365986, -0.9650386511798188, 2.5476426337726172], [1.5424767514304085, -1.2559555411150252, 1.316369312245938], [3.2378450948526947, -0.7569464994859721, 2.5665189880288133], [1.5819802655580877, -2.4596962047296422, 3.6295103703023206], [0.19213413022652728, -2.48231156325931, 3.929922210069819], [2.2889562397950587, -3.5480271463578146, 3.0628733568674056], [2.36862059055278, -1.9965314566037957, 4.912950296825869], [1.746298107969314, 1.383924747436088, 2.5340254928450245], [1.2535978601361861, 0.3461183998849051, 4.378125052116877], [-1.0117971004258965, 1.1239736000383767, 3.432905111838007], [-0.8215282967123819, -0.9363359387499615, 2.788444311143952], [-2.011954187959715, 0.6059455318059196, 1.047594832227991], [-0.40994706586377727, 2.0654783699937855, 1.0587059160250727], [1.5770424436171686, 0.0, 0.0], [2.292718146893918, 1.3915527243580563, 0.0], [1.600521547008252, 2.440721904563901, -0.9357086002340272], [1.3760692761371658, 1.9145350254105713, -2.1257143081021437], [0.4611709133627564, 2.8519467789832516, -0.4188289776134484], [2.404864907693498, 3.4780800111829104, -1.0689303403306858], [3.540059223330468, 1.2295174124846084, -0.4315210556927449], [2.3088468039522523, 1.8960947387583753, 1.2310220414904673], [1.9974224573334847, -0.6906780683055254, 1.053572223549299], [1.9277183224308954, -0.6529932317206267, -1.1102241252095324], [-0.3501493572534735, -1.2838136616209448, 0.08241309473864851], [-0.42668432219275504, 0.4915335256355461, -1.1586058166012287], [2.7837710225593413, -2.7560988038109113, 5.345003477068191]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0301', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
