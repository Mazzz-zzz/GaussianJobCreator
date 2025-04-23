import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0079'
logfile = 'conf/5009017845242299296281_0079.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863813, 0.6217394783082114, -1.2501828803165056], [-0.34660204151390006, 2.1278181305643202, -1.5663863369811297], [-0.7363023803695163, 3.202618946498003, -0.49190218340504926], [-2.2452668735535446, 3.6152156038456105, -0.6031465794771637], [-2.6409584449074375, 4.178902491939337, 0.5247661900346001], [-2.4192802354782867, 4.4675883278057835, -1.6091685125747337], [-3.3507490824222406, 2.134551051026026, -0.9305681716383414], [-3.227636012618887, 1.7823867634696937, -2.3028853705481005], [-3.154695608522511, 1.230044891127082, 0.1412473039565591], [-4.743136573073806, 2.8375865479395794, -0.7125296526432875], [-0.5266683954723087, 2.695152934634473, 0.7195116152968826], [0.0041259347505317985, 4.291386672281077, -0.6555498260725701], [0.980927199419512, 2.156677583916374, -1.7135621949231812], [-0.9145558447621641, 2.4557522019201037, -2.720590366825876], [-2.011954187959717, 0.6042709716797783, -1.0485616399676008], [-0.40994706586377344, -0.11587296658230789, -2.318109697394431], [1.5770424436171668, 0.0, 0.0], [2.292718146893914, 1.3915527243580565, 0.0], [3.7823355744197134, 1.3186147352454554, -0.48079377647244304], [4.419592497958772, 0.34925073248438954, 0.14976318527135613], [3.8473329977166384, 1.1104432852891954, -1.779574085222879], [4.3680104110515785, 2.468434085335201, -0.2056055442087163], [2.293896732220212, 1.8704189044736075, 1.240568989312634], [1.6494649440008828, 2.2352148943143355, -0.8029305726284905], [1.9974224573334802, -0.6906780683055248, 1.0535722235493044], [1.9277183224308945, -0.6529932317206304, -1.1102241252095257], [-0.35014935725347773, 0.7132786644586377, 1.0706086973199298], [-0.4266843221927574, -1.2491488329668494, 0.1536223882885022], [-4.6430901167298995, 3.6110426050797564, -0.14006234271004747]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0079', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
