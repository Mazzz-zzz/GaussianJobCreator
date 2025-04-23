import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0046'
logfile = 'conf/5009017845242299296281_0046.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586382, 0.7718203945763844, 1.1635336229088464], [-0.34660204151390284, 0.2926212946843425, 2.625937724192379], [1.1624535490467065, 0.30841610108309514, 3.054984450098302], [1.9243606525365986, -0.9650386511798188, 2.5476426337726172], [1.5424767514304085, -1.2559555411150252, 1.316369312245938], [3.2378450948526947, -0.7569464994859721, 2.5665189880288133], [1.5819802655580877, -2.4596962047296422, 3.6295103703023206], [0.19213413022652734, -2.4823115632593105, 3.9299222100698197], [2.288956239795058, -3.5480271463578106, 3.062873356867408], [2.36862059055278, -1.9965314566037957, 4.912950296825869], [1.746298107969314, 1.383924747436088, 2.5340254928450245], [1.2535978601361861, 0.3461183998849051, 4.378125052116877], [-1.0117971004258965, 1.1239736000383767, 3.432905111838007], [-0.8215282967123819, -0.9363359387499615, 2.788444311143952], [-2.011954187959715, 0.6059455318059196, 1.047594832227991], [-0.40994706586377727, 2.0654783699937855, 1.0587059160250727], [1.5770424436171686, 0.0, 0.0], [2.292718146893918, 1.3915527243580563, 0.0], [3.782335574419717, 1.3186147352454642, -0.48079377647244653], [4.419592497958775, 0.34925073248439176, 0.1497631852713524], [3.8473329977166393, 1.1104432852892008, -1.7795740852228767], [4.368010411051581, 2.4684340853352067, -0.20560554420871874], [2.2938967322202175, 1.8704189044736075, 1.2405689893126328], [1.6494649440008855, 2.2352148943143364, -0.8029305726284874], [1.9974224573334842, -0.6906780683055251, 1.0535722235492992], [1.9277183224308958, -0.6529932317206266, -1.110224125209532], [-0.3501493572534735, -1.2838136616209448, 0.08241309473864851], [-0.42668432219275504, 0.4915335256355461, -1.1586058166012287], [1.7667412333238892, -1.5746035892577301, 5.541981695719864]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0046', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
