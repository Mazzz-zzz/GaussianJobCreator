import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0470'
logfile = 'conf/5009017845242299296281_0470.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863794, 0.7718203945763795, 1.1635336229088526], [-0.3466020415138972, 0.2926212946843311, 2.6259377241923825], [-0.7363023803695149, -1.1753096862432197, 3.0195004580111413], [-2.2452668735535437, -1.2852675418899158, 3.4324418428267727], [-2.6409584449074375, -2.5439120975868277, 3.3566526229402442], [-2.4192802354782827, -0.8402133530431878, 4.673629241818011], [-3.350749082422239, -0.2613798489209845, 2.313859521682482], [-2.920822010136162, -0.4487371504041443, 0.9712907294388728], [-4.681140306217326, -0.461121463995024, 2.756035109548467], [-2.8983292742780984, 1.1687263611277552, 2.794125255404827], [-0.5266683954723087, -1.9706918044823267, 1.9743151008291795], [0.0041259347505327604, -1.5779705333152534, 4.04422478869365], [0.9809271994195139, 0.40564959980988763, 2.724518672905619], [-0.9145558447621589, 1.1282242700023544, 3.487038975675329], [-2.0119541879597147, 0.6059455318059195, 1.0475948322280018], [-0.4099470658637695, 2.065478369993782, 1.0587059160250882], [1.5770424436171673, 0.0, 0.0], [2.2927181468939173, 1.3915527243580477, 0.0], [1.6005215470082588, 2.440721904563897, -0.9357086002340231], [1.3760692761371665, 1.9145350254105773, -2.1257143081021415], [0.46117091336276195, 2.851946778983256, -0.4188289776134373], [2.404864907693511, 3.478080011182905, -1.0689303403306853], [3.540059223330469, 1.229517412484595, -0.43152105569275456], [2.3088468039522647, 1.8960947387583689, 1.2310220414904665], [1.9974224573334827, -0.690678068305537, 1.0535722235492893], [1.927718322430891, -0.6529932317206271, -1.1102241252095402], [-0.3501493572534775, -1.2838136616209428, 0.08241309473864379], [-0.4266843221927582, 0.49153352563555247, -1.1586058166012252], [-3.64958008471437, 1.778280405539859, 2.7887892029604697]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0470', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
