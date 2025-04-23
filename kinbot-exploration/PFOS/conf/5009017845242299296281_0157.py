import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0157'
logfile = 'conf/5009017845242299296281_0157.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863851, 0.6217394783082126, -1.2501828803165012], [-0.3976197158559567, -0.07566485901595822, -2.6339101198206682], [-1.1233882121466328, 0.4940238789133658, -3.902883790150538], [-0.3678111583333635, 0.11496132015469422, -5.223868290930238], [0.6688104135818955, 0.9167943976383535, -5.394107462597693], [0.05149838568730977, -1.1461400629926153, -5.171879366435432], [-1.475740760567105, 0.27151298615619346, -6.730403024509696], [-2.3402363179501737, -0.8571210633236112, -6.765950477720138], [-1.9223669601590467, 1.6146187425022056, -6.7745475231053245], [-0.3736644517404418, 0.0987422208937877, -7.842212676384806], [-2.351535285418464, -0.013330338437770947, -3.9597873669787527], [-1.1885301103712402, 1.8171105205812559, -3.8278633047836723], [0.9186702629217423, 0.05926285465343824, -2.8183083047232302], [-0.6863393184328396, -1.3666682794365135, -2.523228685784664], [-0.2560445759534336, 1.8735740976390425, -1.391670165756194], [-2.007602477187449, 0.642713061694672, -1.0543092166280623], [1.5770424436171657, 0.0, 0.0], [2.29271814689391, 1.391552724358063, 0.0], [1.6005215470082383, 2.4407219045639015, -0.93570860023403], [1.3760692761371498, 1.9145350254105642, -2.1257143081021495], [0.46117091336273663, 2.8519467789832476, -0.4188289776134527], [2.4048649076934767, 3.4780800111829153, -1.0689303403306913], [3.5400592233304597, 1.2295174124846189, -0.4315210556927501], [2.308846803952244, 1.8960947387583855, 1.2310220414904678], [1.9974224573334896, -0.6906780683055194, 1.0535722235492988], [1.9277183224308958, -0.6529932317206204, -1.1102241252095315], [-0.35014935725347396, 0.7132786644586335, 1.0706086973199336], [-0.4266843221927512, -1.2491488329668539, 0.15362238828850275], [-0.33463778567007546, -0.8220375583199543, -8.136498089245395]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0157', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
