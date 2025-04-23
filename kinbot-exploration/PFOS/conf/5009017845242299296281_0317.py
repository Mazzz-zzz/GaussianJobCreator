import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0317'
logfile = 'conf/5009017845242299296281_0317.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863818, -1.3935598728845981, 0.08664925740765331], [-0.3976197158559555, -2.2432006455416342, 1.382482749991922], [-0.7491833517666546, -1.597752180012834, 2.7686207779580294], [-2.0819874597153203, -0.7736518441984153, 2.7047693461203175], [-1.854244528174576, 0.4066120608580376, 2.1555770167573547], [-2.9956443258697925, -1.4254108214916992, 1.9910877595570866], [-2.8081615356386624, -0.48876666843532235, 4.411571906231566], [-3.7542741288270673, 0.5691492435883607, 4.3214836611193475], [-3.0942747204253434, -1.7648102761531148, 4.954903046827963], [-1.5114899853733603, 0.052549989597533056, 5.1230223938114925], [-0.895313736421604, -2.5657760365450364, 3.6689358763380975], [0.22590414679516813, -0.7802008545701403, 3.1444452075848726], [-1.130548978421929, -3.3520534153271284, 1.2479921638695748], [0.8878652548597824, -2.5746029018458136, 1.3820892141482377], [-0.2560445759534272, -2.1420087660532907, -0.9267276815498183], [-2.007602477187446, -1.2344150958913118, -0.02945123045762922], [1.5770424436171664, 0.0, 0.0], [2.2927181468939173, 1.3915527243580583, 0.0], [2.3410798567223274, 2.059852692894922, 1.4165023767064788], [1.1453916410070415, 2.0331725558915923, 1.975951122830796], [3.202217992284647, 1.442225996159819, 2.1984030628363307], [2.721555559149498, 3.315205977928801, 1.2745358845394033], [1.6292323391939763, 2.2122558673105805, -0.8090479336198823], [3.545586830094381, 1.2600392214310758, -0.42809146886198063], [1.9974224573334851, -0.6906780683055254, 1.0535722235492984], [1.9277183224308962, -0.6529932317206256, -1.1102241252095328], [-0.35014935725347385, 0.5705349971623097, -1.1530217920585801], [-0.4266843221927562, 0.757615307331305, 1.0049834283127268], [-0.7200196046214978, -0.26684471675956484, 4.6674743231981175]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0317', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
