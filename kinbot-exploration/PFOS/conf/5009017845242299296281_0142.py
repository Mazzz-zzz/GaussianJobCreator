import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0142'
logfile = 'conf/5009017845242299296281_0142.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586383, 0.6217394783082169, -1.2501828803165014], [-0.39761971585595723, -0.07566485901595296, -2.6339101198206682], [-0.7491833517666586, -1.5988198371506648, -2.768004365822098], [0.37513150925003946, -2.5102398158058126, -2.164005297273762], [0.8222268376054809, -1.9844600407113075, -1.036990062886352], [-0.09608992658047184, -3.7311518678641797, -1.9268889542995347], [1.8294372022557648, -2.6927750788515747, -3.3359008387581297], [1.4673342643132088, -3.61799689111334, -4.353386498734923], [2.3207428086512274, -1.3887225174009494, -3.5877390601353687], [2.8319974601065314, -3.397143239316472, -2.3463145344316665], [-1.8820373528306378, -1.837820181460768, -2.113395933093811], [-0.8963077255402113, -1.9155699118782823, -4.048025307170632], [-1.13054897842193, 0.595233790028608, -3.526959494450468], [0.887865254859778, 0.09037708117406704, -2.920716124729729], [-0.2560445759534286, 1.8735740976390458, -1.3916701657561938], [-2.007602477187445, 0.6427130616946783, -1.0543092166280612], [1.5770424436171648, 0.0, 0.0], [2.2927181468939155, 1.3915527243580565, 0.0], [2.3410798567223257, 2.0598526928949283, 1.4165023767064775], [1.145391641007037, 2.0331725558915927, 1.975951122830796], [3.202217992284647, 1.4422259961598294, 2.1984030628363236], [2.7215555591494978, 3.3152059779287972, 1.2745358845394015], [1.6292323391939723, 2.2122558673105805, -0.8090479336198859], [3.5455868300943796, 1.2600392214310778, -0.42809146886198063], [1.9974224573334831, -0.690678068305528, 1.0535722235492984], [1.9277183224308954, -0.6529932317206251, -1.1102241252095335], [-0.3501493572534754, 0.7132786644586334, 1.0706086973199334], [-0.42668432219275865, -1.2491488329668496, 0.15362238828850108], [2.6044076003362244, -3.1853916688945003, -1.4301668853872767]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0142', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
