import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0284'
logfile = 'conf/5009017845242299296281_0284.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863816, 0.6217394783082109, -1.2501828803165043], [-0.346602041513903, 2.1278181305643193, -1.566386336981129], [-0.7363023803695204, 3.2026189464980006, -0.4919021834050475], [-0.5090926255515552, 2.6747258360153707, 0.9673064206754453], [-0.4446952607212926, 3.692359865183755, 1.8080060728452527], [-1.5014279159495583, 1.8653940792490031, 1.3266847999265914], [1.0867274347625946, 1.6984177462782482, 1.115104188343934], [1.4232745804972633, 1.6022652107412088, 2.4934806668205645], [0.9696038499690479, 0.585315176193469, 0.24755505082849996], [2.052951777282571, 2.7368576630894763, 0.4305636936294607], [0.01837074119013288, 4.283045989816822, -0.6711713557876985], [-2.0169579203546344, 3.5244552588717037, -0.6214111693177129], [0.980927199419507, 2.156677583916374, -1.7135621949231812], [-0.9145558447621683, 2.4557522019201046, -2.7205903668258746], [-2.011954187959715, 0.6042709716797756, -1.0485616399675983], [-0.4099470658637734, -0.11587296658231044, -2.3181096973944304], [1.577042443617166, 0.0, 0.0], [2.2927181468939137, 1.3915527243580572, 0.0], [2.3410798567223168, 2.0598526928949297, 1.4165023767064777], [1.1453916410070324, 2.0331725558915923, 1.9759511228307929], [3.2022179922846368, 1.4422259961598298, 2.1984030628363285], [2.7215555591494924, 3.3152059779287972, 1.2745358845394077], [1.629232339193974, 2.212255867310578, -0.8090479336198875], [3.5455868300943783, 1.2600392214310787, -0.4280914688619786], [1.997422457333481, -0.6906780683055238, 1.0535722235493046], [1.9277183224308951, -0.6529932317206273, -1.1102241252095277], [-0.35014935725347773, 0.7132786644586361, 1.0706086973199311], [-0.4266843221927561, -1.249148832966851, 0.15362238828850014], [2.733868353261498, 2.276567475100257, -0.07977731598468547]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0284', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
