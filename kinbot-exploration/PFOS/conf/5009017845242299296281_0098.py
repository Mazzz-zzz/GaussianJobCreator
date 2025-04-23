import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0098'
logfile = 'conf/5009017845242299296281_0098.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863837, 0.7718203945763852, 1.1635336229088478], [-0.39761971585595557, 2.3188655045575928, 1.251427369828747], [-0.7491833517666567, 3.1965720171635055, -0.0006164121359324551], [-2.0819874597153216, 2.729224887216827, -0.682382522299651], [-2.5599301491966755, 3.6973095519828596, -1.4445405180545765], [-1.8676581621436081, 1.6463331421236476, -1.4242103063334348], [-3.3945978505277368, 2.2859542163089146, 0.5834090898384298], [-3.3770290100727016, 3.273286151028224, 1.606784873379593], [-4.55380281943126, 1.9133021841908455, -0.13959234492892822], [-2.7366583652515954, 0.968317985308213, 1.1415246928077718], [-0.8953137364216052, 4.4602796920374335, 0.3875592899003024], [0.22590414679516588, 3.113269857861802, -0.8965488436803655], [-1.1305489784219311, 2.7568196252985295, 2.2789673305808957], [0.8878652548597811, 2.4842258206717522, 1.5386269105814903], [-0.25604457595342994, 0.268434668414245, 2.3183978473060143], [-2.007602477187448, 0.5917020341966387, 1.083760447085695], [1.5770424436171646, 0.0, 0.0], [2.292718146893917, 1.3915527243580539, 0.0], [2.341079856722327, 2.05985269289493, 1.4165023767064702], [1.1453916410070402, 2.0331725558915945, 1.9759511228307904], [3.2022179922846457, 1.442225996159828, 2.1984030628363227], [2.721555559149509, 3.3152059779287946, 1.2745358845394004], [1.6292323391939787, 2.2122558673105743, -0.8090479336198904], [3.545586830094382, 1.2600392214310692, -0.42809146886198385], [1.997422457333482, -0.6906780683055271, 1.0535722235493006], [1.9277183224308945, -0.6529932317206282, -1.1102241252095293], [-0.35014935725347623, -1.283813661620943, 0.08241309473865067], [-0.4266843221927585, 0.4915335256355441, -1.158605816601227], [-2.143714196003828, 0.5830362475338164, 0.4812689988613904]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0098', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
