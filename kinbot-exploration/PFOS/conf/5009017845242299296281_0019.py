import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0019'
logfile = 'conf/5009017845242299296281_0019.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586383, -1.3935598728845975, 0.08664925740765078], [-0.39761971585595646, -2.2432006455416356, 1.3824827499919172], [1.087941409756326, -2.6761952078648292, 1.6417185410974606], [1.4730420899455359, -3.9541613281908337, 0.8182592400881729], [0.9614495666061597, -3.8758183253011813, -0.3977331147555902], [2.7953494776680783, -4.070811596542063, 0.7358218224855981], [0.8248630198382829, -5.522044717330672, 1.6202997018164225], [-0.5063280479077278, -5.280040110689663, 2.058243898407448], [1.20042992580093, -6.59887345370361, 0.7808081036958109], [1.7738962514906857, -5.543789448625514, 2.8770437973622625], [1.8947536492841413, -1.6828188678150033, 1.2794354892121673], [1.259810370563856, -2.945657386689846, 2.929461138040973], [-0.7843720924895171, -1.4713296274432142, 2.4020244202192766], [-1.158198376361054, -3.3308565421497387, 1.3537778366861801], [-0.25604457595342883, -2.142008766053284, -0.9267276815498242], [-2.0076024771874454, -1.234415095891309, -0.02945123045763683], [1.5770424436171646, 0.0, 0.0], [2.2927181468939164, 1.3915527243580559, 0.0], [2.3410798567223186, 2.059852692894934, 1.416502376706474], [1.1453916410070355, 2.0331725558915976, 1.9759511228307947], [3.2022179922846465, 1.4422259961598303, 2.1984030628363245], [2.7215555591495, 3.315205977928798, 1.2745358845394037], [1.6292323391939727, 2.2122558673105774, -0.80904793361988], [3.5455868300943765, 1.2600392214310816, -0.4280914688619818], [1.997422457333484, -0.6906780683055264, 1.0535722235492986], [1.9277183224308923, -0.6529932317206255, -1.1102241252095357], [-0.35014935725347696, 0.5705349971623107, -1.1530217920585786], [-0.42668432219275854, 0.7576153073313028, 1.004983428312726], [2.006893193138687, -6.454298323062723, 3.106496192097574]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0019', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
