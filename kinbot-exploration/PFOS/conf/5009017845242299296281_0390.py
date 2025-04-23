import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0390'
logfile = 'conf/5009017845242299296281_0390.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863818, 0.6217394783082115, -1.250182880316503], [-0.3976197158559555, -0.07566485901596232, -2.6339101198206682], [-0.7491833517666546, -1.5988198371506748, -2.768004365822095], [-2.0819874597153203, -1.9555730430184162, -2.0223868238206633], [-1.854244528174576, -2.0700904867547636, -0.72565213419047], [-2.9956443258697925, -1.0116271701948336, -2.229985862019599], [-2.8081615356386624, -3.5761500072006194, -2.629070304503853], [-3.754274128827066, -4.02708925436292, -1.6678441270674467], [-3.0942747204253442, -3.408666773765383, -4.0058220554224055], [-1.5114899853733603, -4.462942531996091, -2.516001570945667], [-0.895313736421604, -1.894503655492401, -4.056495166238398], [0.22590414679516813, -2.3330690032916674, -2.247896363904501], [-1.130548978421929, 0.5952337900285949, -3.5269594944504714], [0.8878652548597823, 0.0903770811740579, -2.920716124729727], [-0.2560445759534272, 1.8735740976390407, -1.3916701657562032], [-2.007602477187446, 0.6427130616946712, -1.0543092166280663], [1.5770424436171664, 0.0, 0.0], [2.2927181468939186, 1.3915527243580539, 0.0], [1.6005215470082594, 2.4407219045639015, -0.9357086002340227], [1.3760692761371767, 1.9145350254105786, -2.125714308102142], [0.46117091336276195, 2.8519467789832538, -0.4188289776134423], [2.4048649076935122, 3.478080011182911, -1.0689303403306736], [3.540059223330471, 1.2295174124846056, -0.43152105569273863], [2.3088468039522505, 1.896094738758371, 1.2310220414904736], [1.9974224573334851, -0.6906780683055258, 1.0535722235493001], [1.9277183224308987, -0.6529932317206291, -1.1102241252095306], [-0.3501493572534738, 0.7132786644586362, 1.0706086973199314], [-0.4266843221927562, -1.249148832966851, 0.15362238828850353], [-1.4815643388024826, -5.107900539286829, -3.236484749933821]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0390', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
