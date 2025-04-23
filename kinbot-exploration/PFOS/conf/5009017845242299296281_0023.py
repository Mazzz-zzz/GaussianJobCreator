import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0023'
logfile = 'conf/5009017845242299296281_0023.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863854, 0.7718203945763855, 1.1635336229088482], [-0.3976197158559566, 2.3188655045575928, 1.2514273698287457], [-0.749183351766658, 3.196572017163504, -0.0006164121359324548], [-2.081987459715323, 2.729224887216826, -0.6823825222996521], [-2.9765439628333943, 2.4254784423769804, 0.24174339068694908], [-2.563171468208477, 3.690357701799602, -1.4656996376850273], [-1.8210733507287753, 1.2133943698125067, -1.7574848436405082], [-1.2114586234507911, 1.6297813413888307, -2.9729776580977867], [-1.303181128745222, 0.19774151930154935, -0.9174885904905876], [-3.332388801547781, 0.872868186472624, -2.041227960534511], [-0.8953137364216069, 4.460279692037433, 0.38755928990030236], [0.22590414679516435, 3.1132698578618005, -0.8965488436803667], [-1.1305489784219311, 2.7568196252985273, 2.2789673305808944], [0.8878652548597796, 2.4842258206717522, 1.538626910581489], [-0.25604457595343105, 0.26843466841424507, 2.3183978473060147], [-2.00760247718745, 0.5917020341966375, 1.0837604470856952], [1.5770424436171646, 0.0, 0.0], [2.2927181468939146, 1.3915527243580543, 0.0], [3.7823355744197134, 1.3186147352454627, -0.4807937764724436], [4.419592497958773, 0.34925073248439165, 0.1497631852713508], [3.847332997716637, 1.110443285289198, -1.779574085222876], [4.368010411051577, 2.4684340853352023, -0.2056055442087178], [2.2938967322202153, 1.870418904473606, 1.2405689893126322], [1.6494649440008793, 2.2352148943143364, -0.8029305726284907], [1.99742245733348, -0.6906780683055269, 1.0535722235493015], [1.9277183224308954, -0.6529932317206273, -1.1102241252095317], [-0.3501493572534779, -1.2838136616209448, 0.08241309473864966], [-0.4266843221927594, 0.49153352563554364, -1.1586058166012267], [-3.4594994221324264, -0.08543365108535307, -2.079475378272115]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0023', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
