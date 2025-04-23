import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0492'
logfile = 'conf/5009017845242299296281_0492.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863848, 0.6217394783082163, -1.2501828803164987], [-0.34660204151390644, 2.127818130564328, -1.5663863369811173], [-0.9873929842445293, 2.7656987490870772, -2.8486010727213635], [-0.18477695124075286, 2.3967981845477464, -4.144481096508231], [0.18832441293677163, 1.1299271526239643, -4.0961581868459], [-0.933236841251377, 2.5972052657952314, -5.225473619786793], [1.3577185656456119, 3.448935518528457, -4.332571345974898], [0.96592075388888, 4.73254537804157, -4.802966085606788], [2.1453375683054725, 3.2344698323419427, -3.1754120882468864], [2.0283879731241963, 2.6593723761190002, -5.518882949150825], [-2.2305068985543404, 2.311247512306057, -2.9794957706412704], [-1.0036314587762365, 4.0874032189906675, -2.733434203914625], [-0.778453490345142, 2.8157771759892394, -0.505705345965246], [0.9731582247379394, 2.2459721722410175, -1.6474522505842994], [-2.0119541879597187, 0.6042709716797802, -1.04856163996759], [-0.40994706586377755, -0.1158729665822974, -2.3181096973944277], [1.5770424436171642, 0.0, 0.0], [2.292718146893915, 1.3915527243580574, 0.0], [1.600521547008249, 2.440721904563901, -0.9357086002340271], [1.3760692761371631, 1.9145350254105762, -2.1257143081021406], [0.46117091336275795, 2.8519467789832573, -0.4188289776134454], [2.404864907693497, 3.47808001118291, -1.0689303403306858], [3.5400592233304664, 1.229517412484604, -0.4315210556927503], [2.308846803952257, 1.8960947387583764, 1.2310220414904702], [1.9974224573334856, -0.6906780683055287, 1.053572223549294], [1.9277183224308947, -0.6529932317206204, -1.110224125209533], [-0.3501493572534748, 0.7132786644586312, 1.0706086973199362], [-0.4266843221927567, -1.2491488329668534, 0.15362238828849817], [1.8085612480262379, 3.073660674333007, -6.365053028739821]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0492', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
