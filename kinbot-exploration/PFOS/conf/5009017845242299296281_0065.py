import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0065'
logfile = 'conf/5009017845242299296281_0065.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863796, -1.3935598728846017, 0.08664925740764856], [-2.2709622836291894, -1.3932545648232224, 0.05367636867327392], [-2.969991788512711, -0.7382819684025497, -1.1888806696804048], [-4.434530032795166, -1.266641279561201, -1.3774259648962097], [-5.118760473256556, -0.4310367516186233, -2.1390661163489013], [-4.423212373540425, -2.4713913694372445, -1.94076349434519], [-5.334436835556408, -1.426961554161535, 0.2615153052569085], [-6.728895866558153, -1.5147133398592965, -0.0034651520876527193], [-4.613374145520475, -2.3678690686870736, 1.0361884090855698], [-5.018585672357864, 0.005273645317363885, 0.8355214079150309], [-3.0171454311974566, 0.5789431357488269, -1.0108747386529298], [-2.2839557356380076, -1.017061933834243, -2.2898193902816013], [-2.6212054717929196, -2.6819879407832867, 0.09005605494849672], [-2.713377980741125, -0.7960712600074846, 1.1536489463716046], [-0.3710451618282775, -1.9490799840121185, 1.2553873021032798], [-0.2455253200204855, -2.1655580448419784, -0.8973245689258842], [1.5770424436171655, 0.0, 0.0], [2.2927181468939155, 1.391552724358058, 0.0], [1.600521547008265, 2.440721904563903, -0.9357086002340264], [1.3760692761371671, 1.9145350254105802, -2.1257143081021397], [0.46117091336276816, 2.8519467789832564, -0.41882897761343973], [2.4048649076934967, 3.478080011182904, -1.0689303403306916], [3.5400592233304655, 1.2295174124846056, -0.4315210556927522], [2.3088468039522647, 1.89609473875837, 1.2310220414904687], [1.9974224573334856, -0.6906780683055247, 1.0535722235492937], [1.9277183224308947, -0.6529932317206142, -1.1102241252095355], [-0.3501493572534782, 0.5705349971623123, -1.1530217920585726], [-0.42668432219275826, 0.7576153073312992, 1.0049834283127321], [-4.214173792050014, 0.3594525689140303, 0.4312427866703521]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0065', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
