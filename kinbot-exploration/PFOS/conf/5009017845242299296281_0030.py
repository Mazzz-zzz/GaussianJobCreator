import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0030'
logfile = 'conf/5009017845242299296281_0030.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863814, -1.3935598728846004, 0.0866492574076522], [-2.2709622836291916, -1.39325456482322, 0.0536763686732763], [-2.9699917885127123, -0.7382819684025475, -1.1888806696804013], [-4.434530032795169, -1.2666412795611983, -1.3774259648962068], [-4.412149490303208, -2.4560673876137575, -1.9529835970041551], [-5.049857983877613, -1.3550224808038633, -0.20164523855451325], [-5.456246996148601, -0.11973325061456042, -2.4555538566829767], [-4.653927843774284, 0.2874984910094085, -3.55686094108644], [-6.738025277307821, -0.7081885753786035, -2.5825244845989586], [-5.582487214586101, 1.0855915838672818, -1.4496721139123991], [-3.0171454311974584, 0.5789431357488257, -1.0108747386529258], [-2.2839557356380094, -1.0170619338342415, -2.2898193902815986], [-2.62120547179292, -2.681987940783285, 0.09005605494849667], [-2.713377980741126, -0.7960712600074821, 1.1536489463716053], [-0.3710451618282789, -1.9490799840121176, 1.2553873021032844], [-0.24552532002048746, -2.165558044841979, -0.897324568925881], [1.5770424436171646, 0.0, 0.0], [2.2927181468939155, 1.3915527243580554, 0.0], [1.6005215470082639, 2.4407219045639046, -0.935708600234023], [1.3760692761371704, 1.9145350254105789, -2.1257143081021406], [0.4611709133627655, 2.8519467789832578, -0.41882897761343807], [2.404864907693504, 3.4780800111829016, -1.0689303403306873], [3.5400592233304695, 1.2295174124846002, -0.43152105569274835], [2.3088468039522643, 1.896094738758368, 1.2310220414904702], [1.9974224573334847, -0.690678068305522, 1.0535722235492966], [1.9277183224308936, -0.6529932317206175, -1.1102241252095344], [-0.3501493572534774, 0.570534997162311, -1.1530217920585715], [-0.426684322192758, 0.7576153073312977, 1.0049834283127321], [-5.478610072590622, 0.7756938829171823, -0.5391034229215996]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0030', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
