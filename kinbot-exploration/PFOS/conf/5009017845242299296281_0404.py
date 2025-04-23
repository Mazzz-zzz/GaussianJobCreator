import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0404'
logfile = 'conf/5009017845242299296281_0404.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863806, 0.621739478308213, -1.2501828803165045], [-0.3466020415139019, 2.1278181305643247, -1.5663863369811237], [-0.9873929842445216, 2.7656987490870715, -2.848601072721373], [-2.4621724959026183, 2.2782191764473745, -3.0653364349826573], [-3.0952492719626608, 2.2335514544131954, -1.906146833019849], [-3.1051333564033308, 3.1002723740922473, -3.8897937793365225], [-2.5252703668709384, 0.5643264923915278, -3.827080392841268], [-3.830277189616423, 0.03608931022527901, -3.6261981393324274], [-1.890670372779976, 0.6520631246236325, -5.089990485647692], [-1.5307227375406445, -0.1587800663564397, -2.842938864656337], [-0.9969507321948188, 4.088639050760528, -2.7117465508075567], [-0.27888546472280534, 2.43143060360588, -3.9194473283466844], [-0.7784534903451408, 2.8157771759892385, -0.5057053459652575], [0.9731582247379454, 2.2459721722410144, -1.647452250584303], [-2.011954187959716, 0.6042709716797786, -1.048561639967597], [-0.40994706586377194, -0.11587296658230209, -2.3181096973944286], [1.5770424436171664, 0.0, 0.0], [2.2927181468939164, 1.391552724358058, 0.0], [3.782335574419716, 1.3186147352454651, -0.4807937764724418], [4.4195924979587735, 0.3492507324843993, 0.149763185271358], [3.8473329977166397, 1.1104432852892023, -1.7795740852228694], [4.368010411051578, 2.4684340853352076, -0.20560554420871363], [2.293896732220214, 1.8704189044736075, 1.2405689893126362], [1.6494649440008817, 2.23521489431434, -0.8029305726284851], [1.9974224573334842, -0.6906780683055249, 1.0535722235493017], [1.927718322430897, -0.6529932317206235, -1.1102241252095266], [-0.35014935725347696, 0.713278664458633, 1.0706086973199327], [-0.4266843221927559, -1.2491488329668505, 0.1536223882885012], [-1.0073517507187166, -0.8186995898076663, -3.3189130440158334]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0404', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
