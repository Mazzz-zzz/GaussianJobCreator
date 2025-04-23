import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0129'
logfile = 'conf/5009017845242299296281_0129.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586383, -1.3935598728845957, 0.08664925740765315], [-0.3976197158559559, -2.2432006455416347, 1.3824827499919194], [1.0879414097563262, -2.6761952078648266, 1.641718541097466], [1.3412747433710008, -3.008180192338106, 3.1534041404090267], [2.4335785182804317, -3.742064612588835, 3.275010324909891], [1.479810301087866, -1.8878195086267913, 3.856610993472195], [-0.08698699647460126, -3.9648101239617124, 3.9061264713717923], [-1.145899229987119, -3.0518111178500336, 4.166066280494233], [-0.23141045828324522, -5.1522080522836475, 3.148112754382456], [0.5934126739618841, -4.351470547390755, 5.2729083249360675], [1.350706773188926, -3.763379719312868, 0.9220132164307516], [1.9101736356829002, -1.696341012641919, 1.2892277976148878], [-0.7843720924895164, -1.471329627443212, 2.402024420219278], [-1.158198376361054, -3.3308565421497387, 1.3537778366861801], [-0.2560445759534282, -2.142008766053285, -0.9267276815498201], [-2.0076024771874468, -1.2344150958913098, -0.029451230457632458], [1.5770424436171635, 0.0, 0.0], [2.292718146893915, 1.3915527243580583, 0.0], [3.782335574419718, 1.318614735245458, -0.48079377647244365], [4.419592497958775, 0.3492507324843924, 0.14976318527135052], [3.8473329977166357, 1.1104432852891988, -1.7795740852228787], [4.368010411051583, 2.468434085335205, -0.20560554420871513], [2.293896732220221, 1.870418904473603, 1.240568989312635], [1.649464944000883, 2.23521489431434, -0.8029305726284817], [1.9974224573334838, -0.6906780683055249, 1.0535722235492988], [1.9277183224308945, -0.6529932317206248, -1.1102241252095322], [-0.3501493572534746, 0.5705349971623115, -1.1530217920585775], [-0.4266843221927573, 0.7576153073313041, 1.0049834283127295], [0.30845357031669957, -5.232640295872877, 5.552730715580002]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0129', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
