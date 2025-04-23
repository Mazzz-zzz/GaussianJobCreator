import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0396'
logfile = 'conf/5009017845242299296281_0396.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863821, -1.3935598728846001, 0.08664925740765095], [-0.3976197158559556, -2.2432006455416387, 1.3824827499919163], [-0.7491833517666582, -1.5977521800128405, 2.768620777958024], [0.3751315092500395, -0.6189636534602693, 3.255934098715885], [-0.10815375808507581, 0.19767795751880263, 4.175735040589357], [1.3965384052652352, -1.3030408487940415, 3.7635425936703073], [1.0485429431935351, 0.42508470237727286, 1.8494692387133023], [1.7447219915669403, 1.5310634161303949, 2.410340725901663], [1.6092624577838985, -0.46923022843379275, 0.9054611867699176], [-0.32333237119353386, 0.9345024441773525, 1.2672370548521799], [-1.8820373528306364, -0.9113444755835778, 2.6482969312796567], [-0.8963077255402097, -2.5479077952329368, 3.6829448599970274], [-1.1305489784219278, -3.3520534153271333, 1.2479921638695664], [0.8878652548597813, -2.5746029018458185, 1.3820892141482315], [-0.25604457595342656, -2.142008766053287, -0.9267276815498245], [-2.007602477187446, -1.2344150958913118, -0.02945123045763415], [1.5770424436171644, 0.0, 0.0], [2.2927181468939137, 1.391552724358056, 0.0], [2.3410798567223177, 2.059852692894931, 1.4165023767064737], [1.1453916410070324, 2.033172555891597, 1.9759511228307929], [3.202217992284644, 1.4422259961598303, 2.1984030628363267], [2.7215555591494924, 3.3152059779288012, 1.2745358845394028], [1.6292323391939711, 2.212255867310577, -0.8090479336198844], [3.5455868300943747, 1.2600392214310787, -0.428091468861986], [1.9974224573334838, -0.6906780683055278, 1.0535722235492968], [1.9277183224308958, -0.6529932317206237, -1.1102241252095326], [-0.3501493572534777, 0.5705349971623089, -1.153021792058581], [-0.42668432219276087, 0.7576153073313022, 1.0049834283127268], [-1.0346741316976078, 0.3203281612820711, 1.4969003664390144]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0396', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
