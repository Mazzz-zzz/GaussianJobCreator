import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0149'
logfile = 'conf/5009017845242299296281_0149.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863868, -1.3935598728845988, 0.08664925740765335], [-0.34660204151390495, -2.420439425248668, -1.0595513872112574], [1.1624535490467045, -2.7999021924931133, -1.260396046575038], [1.3229555267971613, -4.171696370783135, -2.0033664407024196], [0.43579710458927745, -4.260501353184645, -2.978848655857203], [2.548490142211768, -4.284594543793282, -2.5075677392970843], [1.0607552256216175, -5.626043159278573, -0.8467378261629634], [2.2341818168842686, -5.7852102249496244, -0.059231071948513785], [-0.24309140467177717, -5.494666535316873, -0.3100061250338068], [1.0184448819506249, -6.771214053041932, -1.9271711136475143], [1.7500514572469164, -1.8550381832545386, -1.9890607533025857], [1.763328214900581, -2.894514914439559, -0.08123199086867498], [-1.0117971004259, -3.534969835652362, -0.7430628651027226], [-0.8215282967123799, -1.9466956411138836, -2.2051128650057956], [-2.0119541879597196, -1.210216503485698, 0.0009668077395966854], [-0.4099470658637843, -1.949605403411482, 1.259403781369355], [1.5770424436171646, 0.0, 0.0], [2.292718146893919, 1.3915527243580537, 0.0], [1.60052154700826, 2.4407219045638975, -0.9357086002340302], [1.376069276137173, 1.914535025410578, -2.125714308102146], [0.4611709133627615, 2.8519467789832555, -0.4188289776134505], [2.4048649076935122, 3.4780800111829038, -1.0689303403306814], [3.5400592233304686, 1.2295174124846004, -0.43152105569274085], [2.3088468039522585, 1.896094738758367, 1.231022041490469], [1.9974224573334802, -0.6906780683055256, 1.0535722235493001], [1.927718322430895, -0.6529932317206284, -1.1102241252095275], [-0.35014935725347374, 0.5705349971623063, -1.153021792058585], [-0.42668432219276103, 0.7576153073313048, 1.0049834283127228], [1.8912046102274624, -7.178239586343057, -2.0198159118559063]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0149', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
