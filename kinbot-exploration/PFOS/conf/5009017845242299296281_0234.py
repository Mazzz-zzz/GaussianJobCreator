import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0234'
logfile = 'conf/5009017845242299296281_0234.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863816, 0.621739478308212, -1.2501828803165038], [-0.3466020415139052, 2.1278181305643233, -1.5663863369811226], [-0.9873929842445249, 2.765698749087069, -2.848601072721373], [-2.462172495902622, 2.278219176447373, -3.065336434982655], [-3.095249271962663, 2.2335514544131936, -1.9061468330198459], [-3.1051333564033357, 3.100272374092246, -3.8897937793365203], [-2.525270366870942, 0.5643264923915257, -3.8270803928412653], [-2.239776335865652, 0.6837363734396475, -5.215124140847219], [-1.8335865604435533, -0.30267292708071303, -2.946601948920613], [-4.066604391334008, 0.2932738173956116, -3.6499857790918573], [-0.9969507321948232, 4.088639050760528, -2.7117465508075567], [-0.27888546472281034, 2.4314306036058784, -3.919447328346685], [-0.7784534903451427, 2.8157771759892376, -0.5057053459652561], [0.9731582247379439, 2.245972172241014, -1.647452250584306], [-2.0119541879597174, 0.6042709716797778, -1.0485616399675957], [-0.40994706586377433, -0.11587296658230409, -2.3181096973944273], [1.5770424436171646, 0.0, 0.0], [2.292718146893914, 1.3915527243580577, 0.0], [1.600521547008249, 2.440721904563901, -0.9357086002340301], [1.3760692761371647, 1.9145350254105789, -2.1257143081021432], [0.4611709133627553, 2.85194677898325, -0.4188289776134523], [2.404864907693492, 3.4780800111829135, -1.0689303403306858], [3.5400592233304646, 1.2295174124846102, -0.4315210556927461], [2.3088468039522505, 1.896094738758376, 1.2310220414904687], [1.9974224573334864, -0.690678068305524, 1.0535722235493006], [1.9277183224308967, -0.6529932317206218, -1.1102241252095306], [-0.3501493572534756, 0.7132786644586333, 1.0706086973199334], [-0.4266843221927561, -1.2491488329668508, 0.15362238828850236], [-4.420134739927349, 0.8255298087287122, -2.923568379100991]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0234', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
