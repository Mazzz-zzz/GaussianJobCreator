import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0366'
logfile = 'conf/5009017845242299296281_0366.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586388, 0.6217394783082169, -1.250182880316497], [-0.3466020415139042, 2.1278181305643313, -1.5663863369811144], [-0.9873929842445259, 2.7656987490870835, -2.84860107272136], [-2.4621724959026237, 2.2782191764473905, -3.0653364349826413], [-2.461597984833636, 1.071414784733832, -3.603968752203673], [-3.114143520324482, 2.2457738046910283, -1.9065521057311154], [-3.414181650933439, 3.422830080870188, -4.207793338908967], [-4.567822843373354, 2.731386512418116, -4.66979853212878], [-3.4479075882278027, 4.692976610261773, -3.5825356906741477], [-2.3659963988850503, 3.4966168736067074, -5.380996567857356], [-0.9969507321948222, 4.08863905076054, -2.7117465508075385], [-0.27888546472281345, 2.4314306036058935, -3.919447328346675], [-0.7784534903451381, 2.815777175989243, -0.5057053459652442], [0.9731582247379437, 2.245972172241021, -1.647452250584299], [-2.011954187959716, 0.6042709716797869, -1.0485616399675868], [-0.40994706586377944, -0.11587296658229425, -2.3181096973944264], [1.5770424436171655, 0.0, 0.0], [2.292718146893919, 1.3915527243580592, 0.0], [3.7823355744197174, 1.3186147352454622, -0.4807937764724468], [4.4195924979587735, 0.34925073248439176, 0.14976318527134946], [3.847332997716636, 1.110443285289202, -1.7795740852228763], [4.368010411051582, 2.468434085335203, -0.20560554420871757], [2.2938967322202215, 1.870418904473603, 1.2405689893126357], [1.6494649440008837, 2.2352148943143417, -0.802930572628482], [1.9974224573334887, -0.6906780683055267, 1.0535722235492961], [1.9277183224308925, -0.652993231720621, -1.11022412520953], [-0.3501493572534711, 0.7132786644586333, 1.0706086973199374], [-0.42668432219275587, -1.2491488329668485, 0.15362238828849983], [-1.4772738756603097, 3.3050782264702963, -5.050167896091742]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0366', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
