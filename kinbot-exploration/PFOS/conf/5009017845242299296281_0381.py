import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0381'
logfile = 'conf/5009017845242299296281_0381.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863824, 0.621739478308219, -1.2501828803165025], [-0.3976197158559557, -0.07566485901595067, -2.63391011982067], [-0.7491833517666568, -1.5988198371506621, -2.768004365822102], [0.3751315092500395, -2.5102398158058112, -2.1640052972737647], [0.8222268376054805, -1.984460040711305, -1.036990062886355], [-0.09608992658046996, -3.7311518678641766, -1.926888954299539], [1.8294372022557681, -2.692775078851572, -3.335900838758132], [2.112852701691653, -1.4163542730249474, -3.895236685045862], [2.7882050813003083, -3.501477040451924, -2.678540919791226], [1.128311689997794, -3.580128033894059, -4.432062496109529], [-1.8820373528306351, -1.8378201814607642, -2.113395933093816], [-0.8963077255402097, -1.915569911878278, -4.048025307170637], [-1.1305489784219276, 0.5952337900286115, -3.5269594944504696], [0.8878652548597809, 0.09037708117407094, -2.920716124729729], [-0.25604457595342717, 1.8735740976390478, -1.391670165756193], [-2.0076024771874454, 0.6427130616946796, -1.0543092166280652], [1.5770424436171644, 0.0, 0.0], [2.292718146893914, 1.391552724358058, 0.0], [2.3410798567223163, 2.0598526928949297, 1.416502376706479], [1.1453916410070302, 2.03317255589159, 1.9759511228307955], [3.202217992284642, 1.442225996159832, 2.1984030628363245], [2.7215555591494875, 3.3152059779288, 1.2745358845394028], [1.629232339193965, 2.2122558673105805, -0.8090479336198851], [3.545586830094379, 1.2600392214310816, -0.42809146886198], [1.9974224573334811, -0.6906780683055271, 1.0535722235492988], [1.9277183224308958, -0.6529932317206242, -1.1102241252095344], [-0.35014935725347546, 0.7132786644586324, 1.070608697319934], [-0.4266843221927598, -1.2491488329668488, 0.1536223882884976], [0.7947630307666215, -3.0251147246835326, -5.150858646085492]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0381', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
