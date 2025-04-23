import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0237'
logfile = 'conf/5009017845242299296281_0237.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586385, 0.6217394783082189, -1.2501828803164996], [-0.3466020415139042, 2.127818130564328, -1.5663863369811188], [-0.9873929842445244, 2.765698749087077, -2.8486010727213684], [-1.0274356123761041, 4.331050317691674, -2.76110128071334], [-1.1749864386042272, 4.844696265053522, -3.969777638352199], [-2.0336420090655936, 4.725901990131815, -1.9861555506229136], [0.5505558709515904, 5.028814916256168, -2.023042492641287], [1.6529980492054575, 4.333353033447739, -2.591712010883138], [0.41632863418113125, 6.438545282787703, -2.023380480580311], [0.3487395885195963, 4.531007578835056, -0.5425018235733804], [-0.2614927667443652, 2.421588654926837, -3.9086240299426986], [-2.230877915023897, 2.326974815430626, -2.9958978427152325], [-0.778453490345141, 2.8157771759892407, -0.5057053459652502], [0.9731582247379432, 2.2459721722410166, -1.6474522505842988], [-2.0119541879597183, 0.604270971679784, -1.0485616399675926], [-0.4099470658637806, -0.11587296658229958, -2.31810969739443], [1.5770424436171635, 0.0, 0.0], [2.2927181468939133, 1.3915527243580565, 0.0], [2.3410798567223186, 2.0598526928949283, 1.4165023767064782], [1.1453916410070295, 2.0331725558915914, 1.9759511228307913], [3.202217992284642, 1.4422259961598276, 2.1984030628363267], [2.721555559149495, 3.3152059779287995, 1.2745358845394061], [1.6292323391939658, 2.212255867310578, -0.8090479336198867], [3.5455868300943747, 1.2600392214310794, -0.4280914688619842], [1.997422457333482, -0.6906780683055279, 1.0535722235493004], [1.9277183224308894, -0.6529932317206253, -1.1102241252095282], [-0.35014935725347723, 0.7132786644586345, 1.0706086973199351], [-0.42668432219276314, -1.2491488329668512, 0.1536223882884979], [0.6931982976311326, 5.184885421822911, 0.08180126857238419]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0237', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
