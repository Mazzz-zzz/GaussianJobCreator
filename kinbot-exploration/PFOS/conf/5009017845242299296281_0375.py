import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0375'
logfile = 'conf/5009017845242299296281_0375.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863848, -1.393559872884596, 0.08664925740765565], [-0.346602041513907, -2.4204394252486674, -1.0595513872112532], [-0.9873929842445301, -3.849810268767846, -0.9708648395635582], [-0.1847769512407539, -4.787625007354387, -0.0034475673086427476], [0.1883244129367716, -4.112340624040131, 1.0695334748247955], [-0.9332368412513781, -5.823995534438403, 0.36349107087203164], [1.3577185656456106, -5.476584608587031, -0.8205801020726293], [2.2189737591343435, -5.968557133194488, 0.19855361561609927], [0.9276165596370575, -6.240330490564647, -1.9328126604153857], [1.9474063794356817, -4.12368813154154, -1.3705637695341215], [-2.2305068985543404, -3.735942783996663, -0.5118511747699831], [-1.0036314587762374, -4.410925069658697, -2.1730779211988813], [-0.7784534903451426, -1.845842264430128, -2.1856818928204578], [0.9731582247379391, -2.5497215866483613, -1.1213428320614816], [-2.0119541879597183, -1.210216503485695, 0.0009668077396042073], [-0.4099470658637772, -1.9496054034114785, 1.2594037813693577], [1.5770424436171644, 0.0, 0.0], [2.2927181468939186, 1.3915527243580519, 0.0], [1.6005215470082639, 2.440721904563897, -0.9357086002340244], [1.3760692761371696, 1.9145350254105797, -2.1257143081021415], [0.46117091336276417, 2.8519467789832555, -0.4188289776134432], [2.4048649076935074, 3.4780800111829038, -1.0689303403306802], [3.540059223330467, 1.229517412484598, -0.4315210556927434], [2.308846803952256, 1.8960947387583684, 1.2310220414904733], [1.9974224573334802, -0.690678068305527, 1.0535722235493012], [1.9277183224308918, -0.6529932317206272, -1.1102241252095308], [-0.35014935725347507, 0.5705349971623106, -1.1530217920585806], [-0.42668432219275526, 0.7576153073313066, 1.004983428312729], [1.2419200179155014, -3.471772713899673, -1.485705719516946]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0375', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
