import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0319'
logfile = 'conf/5009017845242299296281_0319.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863808, 0.7718203945763867, 1.1635336229088478], [-2.270962283629191, 0.7431123812655668, 1.1797556627389008], [-2.9970239643018997, 1.4656877041954839, 2.368157397369682], [-4.4759542789311215, 1.8399334351693144, 2.0044399021614114], [-5.036363599701665, 0.8513196024881734, 1.3298392602242735], [-5.178548597123876, 2.076901973905839, 3.1085717675857136], [-4.5636267180482815, 3.3871491523351, 0.9460919446874013], [-3.5410856494537764, 3.3105214880315548, -0.039285633660606965], [-5.931513953723526, 3.6087609985547004, 0.6543223922322354], [-4.133463285573643, 4.435597562461193, 2.0398277365835504], [-3.015076885801341, 0.6515380066656313, 3.4197399735582468], [-2.3510956819310196, 2.582573574372415, 2.6776285729604337], [-2.6427825570536125, 1.3431159399560095, 0.04548369333676643], [-2.667789377892033, -0.523270637705672, 1.1465664052870899], [-0.3710451618282821, 2.061737287215915, 1.0602591291106096], [-0.24552532002048844, 0.3056731502912632, 2.324090564665853], [1.5770424436171666, 0.0, 0.0], [2.292718146893917, 1.3915527243580592, 0.0], [1.6005215470082468, 2.4407219045638984, -0.935708600234028], [1.3760692761371667, 1.9145350254105744, -2.1257143081021463], [0.4611709133627522, 2.8519467789832507, -0.41882897761345295], [2.404864907693493, 3.4780800111829127, -1.0689303403306831], [3.5400592233304664, 1.2295174124846127, -0.4315210556927468], [2.308846803952252, 1.896094738758377, 1.2310220414904687], [1.9974224573334824, -0.6906780683055258, 1.0535722235492968], [1.9277183224308976, -0.6529932317206254, -1.1102241252095322], [-0.350149357253478, -1.283813661620943, 0.08241309473864838], [-0.426684322192755, 0.49153352563554853, -1.1586058166012267], [-3.184449598276523, 4.6137189821968425, 1.979759709729926]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0319', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
