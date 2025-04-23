import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0465'
logfile = 'conf/5009017845242299296281_0465.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586382, -1.3935598728845968, 0.08664925740765447], [-2.270962283629191, -1.3932545648232155, 0.053676368673277366], [-2.9970239643019005, -2.7837283183799246, 0.08524408716294152], [-2.2704853996782592, -3.797387686862249, 1.036149738531485], [-1.8789938252023834, -3.179148869065836, 2.1366119388641853], [-3.086123191320812, -4.797903999452992, 1.3564654291201523], [-0.7559231095728453, -4.5519317793882825, 0.2247142124792245], [-1.1920897005380096, -5.548168294996848, -0.6916386222432571], [0.10363285767357697, -3.4801350358769074, -0.1184651158864186], [-0.14900829291124634, -5.262725983416888, 1.4923768784298173], [-4.239914773111563, -2.609835763869494, 0.5256704367020001], [-3.02305458293385, -3.305441500441076, -1.134414702846077], [-2.6427825570536156, -0.7109480038655792, 1.1404306775613309], [-2.667789377892032, -0.7313203152515892, -1.0264488679511423], [-0.3710451618282834, -1.9490799840121071, 1.255387302103292], [-0.245525320020488, -2.1655580448419838, -0.8973245689258694], [1.5770424436171655, 0.0, 0.0], [2.2927181468939173, 1.3915527243580572, 0.0], [1.6005215470082543, 2.4407219045638993, -0.9357086002340305], [1.3760692761371702, 1.914535025410578, -2.125714308102144], [0.4611709133627524, 2.851946778983251, -0.4188289776134527], [2.4048649076934976, 3.4780800111829127, -1.068930340330687], [3.54005922333046, 1.2295174124846109, -0.43152105569274335], [2.308846803952242, 1.896094738758376, 1.231022041490471], [1.9974224573334833, -0.6906780683055223, 1.053572223549302], [1.927718322430896, -0.6529932317206273, -1.1102241252095306], [-0.35014935725347235, 0.5705349971623055, -1.1530217920585832], [-0.42668432219275565, 0.7576153073313071, 1.0049834283127224], [-0.4506169289251351, -4.824129539949333, 2.3002306877019203]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0465', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
