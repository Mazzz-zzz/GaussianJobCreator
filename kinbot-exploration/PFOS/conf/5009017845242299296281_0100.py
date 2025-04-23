import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0100'
logfile = 'conf/5009017845242299296281_0100.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863868, 0.6217394783082142, -1.2501828803164987], [-0.3466020415139019, 2.1278181305643287, -1.5663863369811182], [-0.9873929842445243, 2.7656987490870772, -2.8486010727213666], [-2.4621724959026214, 2.278219176447384, -3.0653364349826466], [-2.461597984833634, 1.0714147847338245, -3.6039687522036767], [-3.1141435203244794, 2.245773804691025, -1.906552105731123], [-3.414181650933437, 3.42283008087018, -4.207793338908976], [-2.567225413766637, 3.7702224011462997, -5.29610152869137], [-4.71012561367275, 2.871455404216755, -4.355685133924346], [-3.530576733854341, 4.6621295945134875, -3.2428417701946155], [-0.996950732194819, 4.088639050760535, -2.711746550807548], [-0.27888546472281034, 2.4314306036058846, -3.9194473283466795], [-0.7784534903451366, 2.815777175989242, -0.5057053459652504], [0.973158224737945, 2.245972172241018, -1.6474522505843026], [-2.011954187959716, 0.6042709716797844, -1.0485616399675897], [-0.40994706586377694, -0.11587296658229841, -2.3181096973944273], [1.5770424436171655, 0.0, 0.0], [2.2927181468939186, 1.39155272435806, 0.0], [1.6005215470082503, 2.440721904563903, -0.9357086002340288], [1.3760692761371665, 1.9145350254105762, -2.125714308102142], [0.46117091336275906, 2.8519467789832538, -0.4188289776134485], [2.4048649076934976, 3.4780800111829118, -1.0689303403306858], [3.5400592233304664, 1.2295174124846076, -0.43152105569274896], [2.3088468039522554, 1.8960947387583755, 1.231022041490468], [1.9974224573334887, -0.6906780683055248, 1.0535722235492992], [1.9277183224308923, -0.6529932317206233, -1.110224125209528], [-0.35014935725347224, 0.7132786644586361, 1.0706086973199351], [-0.4266843221927561, -1.2491488329668483, 0.15362238828850205], [-4.38084510724131, 5.106479229955815, -3.3675794204696765]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0100', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
