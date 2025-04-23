import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0121'
logfile = 'conf/5009017845242299296281_0121.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863813, -1.3935598728846004, 0.08664925740765345], [-0.34660204151389995, -2.420439425248667, -1.0595513872112567], [-0.7363023803695175, -2.0273092602547975, -2.5275982746061056], [-0.5090926255515525, -0.4996509844589459, -2.8000337324856055], [0.6213430942919618, -0.10684176608510124, -2.239393077217541], [-0.4596848581732949, -0.26494129642363706, -4.108227248529034], [-1.8990683064763114, 0.5480559603545139, -2.098648606223012], [-3.018501756452613, 0.442481731255444, -2.969368906339119], [-1.94865260897711, 0.2875063150960567, -0.7076015002811574], [-1.2380964843927418, 1.9640389989312392, -2.2952469650248015], [0.018370741190135775, -2.7227744393129925, -3.373640954864593], [-2.016957920354632, -2.300385488260377, -2.741562204025706], [0.9809271994195116, -2.562327183726282, -1.0109564779824425], [-0.9145558447621639, -3.583976471922482, -0.7664486088494542], [-2.011954187959714, -1.210216503485702, 0.0009668077396031379], [-0.4099470658637693, -1.9496054034114834, 1.2594037813693562], [1.5770424436171675, 0.0, 0.0], [2.2927181468939155, 1.3915527243580579, 0.0], [3.7823355744197142, 1.318614735245466, -0.48079377647244537], [4.41959249795878, 0.3492507324844002, 0.14976318527134158], [3.847332997716636, 1.1104432852892012, -1.7795740852228819], [4.368010411051581, 2.4684340853352067, -0.20560554420872112], [2.293896732220215, 1.8704189044736084, 1.2405689893126355], [1.6494649440008788, 2.235214894314339, -0.8029305726284843], [1.9974224573334871, -0.6906780683055253, 1.0535722235492986], [1.9277183224308962, -0.6529932317206243, -1.1102241252095322], [-0.35014935725347646, 0.5705349971623088, -1.1530217920585797], [-0.4266843221927541, 0.7576153073313024, 1.004983428312729], [-1.4776687565480118, 2.5511131313918316, -1.5645545015998614]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0121', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
