import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0303'
logfile = 'conf/5009017845242299296281_0303.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863834, 0.6217394783082172, -1.2501828803165005], [-0.3466020415139029, 2.127818130564326, -1.566386336981119], [1.1624535490467056, 2.4914860914100188, -1.7945884035232569], [1.3229555267971613, 3.820814416129082, -2.6111118136223213], [2.537966839189117, 4.31069906016461, -2.437150641590282], [1.1168882370515762, 3.596724293045492, -3.905800938233435], [0.0882698971352805, 5.125787773594594, -2.0686551210406576], [0.5381926157038321, 6.388500124917685, -2.543584666280777], [-1.2025433687641451, 4.598922353138966, -2.3167895865108727], [0.34402262312628656, 5.0496785798237696, -0.5164500561150444], [1.750051457246916, 2.650096233657925, -0.611979815037262], [1.7633282149005842, 1.517606424912047, -2.4661074521032558], [-1.0117971004258974, 2.4109962356139945, -2.68984224673528], [-0.8215282967123809, 2.883031579863846, -0.5833314461381526], [-2.0119541879597174, 0.6042709716797827, -1.04856163996759], [-0.4099470658637775, -0.11587296658229856, -2.3181096973944304], [1.577042443617167, 0.0, 0.0], [2.292718146893916, 1.3915527243580552, 0.0], [1.6005215470082508, 2.4407219045638984, -0.9357086002340284], [1.3760692761371627, 1.914535025410574, -2.1257143081021437], [0.46117091336275684, 2.851946778983252, -0.41882897761344706], [2.404864907693494, 3.4780800111829118, -1.0689303403306871], [3.5400592233304655, 1.229517412484605, -0.4315210556927491], [2.3088468039522527, 1.8960947387583738, 1.2310220414904696], [1.9974224573334851, -0.6906780683055272, 1.0535722235492984], [1.9277183224308951, -0.6529932317206235, -1.1102241252095353], [-0.35014935725347357, 0.7132786644586349, 1.0706086973199356], [-0.4266843221927547, -1.2491488329668514, 0.1536223882885002], [0.9755142919476889, 5.73092422387377, -0.24610762258294247]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0303', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
