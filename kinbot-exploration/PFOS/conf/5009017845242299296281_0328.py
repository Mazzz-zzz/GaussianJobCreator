import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0328'
logfile = 'conf/5009017845242299296281_0328.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863841, -1.393559872884598, 0.08664925740765081], [-0.39761971585595707, -2.243200645541637, 1.3824827499919152], [1.0879414097563251, -2.6761952078648292, 1.6417185410974606], [2.1006320219594934, -1.5474952427548003, 1.241619463720128], [3.2636602341713647, -1.7574463535228209, 1.833048116867345], [2.2719309703990027, -1.525272508610281, -0.07711586260721295], [1.4886616557511099, 0.15148205615249102, 1.7521405110780743], [0.48999533174814114, 0.5612428760438004, 0.826293519618688], [1.3031116920812955, 0.11707668162720598, 3.155616286727196], [2.80710863078977, 0.963856986862653, 1.4652830974139466], [1.2395907752722295, -2.9383364670191545, 2.9368154419869685], [1.3717373297746493, -3.759054226052555, 0.9294643561329469], [-0.7843720924895173, -1.4713296274432146, 2.4020244202192775], [-1.1581983763610537, -3.3308565421497422, 1.3537778366861781], [-0.25604457595342883, -2.142008766053283, -0.9267276815498262], [-2.007602477187449, -1.2344150958913087, -0.02945123045763408], [1.5770424436171655, 0.0, 0.0], [2.2927181468939164, 1.3915527243580568, 0.0], [1.6005215470082472, 2.4407219045638984, -0.9357086002340265], [1.3760692761371598, 1.914535025410576, -2.1257143081021463], [0.4611709133627515, 2.851946778983254, -0.41882897761344484], [2.4048649076934985, 3.4780800111829118, -1.0689303403306867], [3.5400592233304606, 1.229517412484606, -0.43152105569275107], [2.3088468039522536, 1.8960947387583755, 1.2310220414904651], [1.997422457333482, -0.6906780683055266, 1.0535722235492968], [1.9277183224308936, -0.652993231720624, -1.1102241252095342], [-0.3501493572534735, 0.5705349971623116, -1.1530217920585804], [-0.42668432219275926, 0.7576153073313046, 1.0049834283127281], [2.918220908209195, 1.6665566814273292, 2.12089699885395]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0328', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
