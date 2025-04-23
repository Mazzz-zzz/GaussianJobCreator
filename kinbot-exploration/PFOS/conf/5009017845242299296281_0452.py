import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0452'
logfile = 'conf/5009017845242299296281_0452.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863828, 0.6217394783082142, -1.2501828803164985], [-0.397619715855957, -0.07566485901595642, -2.6339101198206665], [1.0879414097563236, -0.08367235852191894, -3.1385123060458504], [1.473042089945533, 1.2684473752977068, -3.833533780919367], [1.000473775901915, 1.2895411745485026, -5.067525751424903], [0.9773040172258548, 2.2971743117836962, -3.1517239144922664], [3.332923574810669, 1.5027949558458802, -3.9206280107056743], [3.791719895996161, 1.9011060831692883, -2.63482663050459], [3.8528498715265944, 0.40081741942365573, -4.642208863314768], [3.3600157795274916, 2.7542140279606304, -4.876535405906718], [1.8947536492841397, -0.2666142022536058, -2.097081634101644], [1.2598103705638541, -1.064159071597835, -4.0157446967391754], [-0.7843720924895167, -1.344549354698879, -2.475221044816154], [-1.158198376361056, 0.49302227342429294, -3.561495300206356], [-0.25604457595342917, 1.8735740976390431, -1.3916701657561907], [-2.007602477187446, 0.6427130616946728, -1.0543092166280583], [1.5770424436171657, 0.0, 0.0], [2.292718146893915, 1.3915527243580559, 0.0], [2.341079856722322, 2.05985269289493, 1.4165023767064724], [1.1453916410070375, 2.033172555891591, 1.9759511228307938], [3.2022179922846457, 1.4422259961598272, 2.198403062836324], [2.7215555591495, 3.315205977928798, 1.2745358845394001], [1.6292323391939738, 2.212255867310576, -0.8090479336198863], [3.5455868300943783, 1.260039221431076, -0.428091468861984], [1.997422457333487, -0.6906780683055249, 1.0535722235492984], [1.9277183224308962, -0.6529932317206225, -1.1102241252095335], [-0.35014935725347246, 0.713278664458632, 1.0706086973199334], [-0.4266843221927523, -1.2491488329668534, 0.1536223882885027], [2.555001092695247, 2.7854138727430664, -5.412203503114244]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0452', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
