import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0136'
logfile = 'conf/5009017845242299296281_0136.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863845, 0.7718203945763872, 1.1635336229088507], [-0.39761971585595557, 2.3188655045575923, 1.2514273698287481], [1.0879414097563276, 2.759867566386744, 1.4967937649483836], [2.1006320219594956, 1.8490216187922401, 0.7193604606011568], [2.2762575910635965, 0.7136703779757695, 1.3725336095673062], [1.6464483143989894, 1.5907192761912547, -0.5037082619759717], [3.775607346559006, 2.6734660435613034, 0.528540169319106], [3.678229033555621, 3.6442164663676326, -0.5061540312000352], [4.248696567169458, 2.9504295588454275, 1.8342333079208424], [4.578856613177445, 1.4332091642747282, -0.016525073443874234], [1.2395907752722297, 4.012525012496718, 1.0762663043113097], [1.371737329774654, 2.6844668573495585, 2.7907042758982876], [-0.7843720924895153, 2.8158789821420935, 0.07319662459687713], [-1.1581983763610517, 2.8378342687254476, 2.2077174635201775], [-0.2560445759534292, 0.2684346684142433, 2.3183978473060174], [-2.007602477187448, 0.5917020341966387, 1.083760447085695], [1.5770424436171642, 0.0, 0.0], [2.292718146893914, 1.3915527243580537, 0.0], [1.6005215470082517, 2.4407219045638993, -0.9357086002340298], [1.3760692761371645, 1.9145350254105715, -2.125714308102145], [0.46117091336275773, 2.8519467789832493, -0.4188289776134525], [2.404864907693498, 3.478080011182908, -1.0689303403306885], [3.5400592233304637, 1.229517412484604, -0.431521055692749], [2.3088468039522514, 1.8960947387583738, 1.2310220414904682], [1.997422457333481, -0.6906780683055286, 1.0535722235493008], [1.9277183224308918, -0.6529932317206286, -1.1102241252095302], [-0.35014935725347895, -1.2838136616209421, 0.08241309473865062], [-0.42668432219276153, 0.4915335256355454, -1.158605816601225], [4.154297978632877, 0.609443973387857, 0.26118568042203105]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0136', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
