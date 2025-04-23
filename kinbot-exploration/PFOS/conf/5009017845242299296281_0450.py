import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0450'
logfile = 'conf/5009017845242299296281_0450.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863837, -1.393559872884597, 0.08664925740765324], [-2.2709622836291934, -1.3932545648232137, 0.053676368673281016], [-2.969991788512716, -0.7382819684025406, -1.1888806696803962], [-3.0567398040787155, 0.8216102792205098, -1.0517817396959934], [-4.057927071959131, 1.1471549786054227, -0.2529409255073739], [-1.9208515428855029, 1.306979118885896, -0.5587553832825151], [-3.34387825427727, 1.6524617020879773, -2.709762666803323], [-2.111023082257976, 1.6689337329433178, -3.4184429930681137], [-4.551026485388037, 1.1257236879582997, -3.230048150211123], [-3.6449190406827388, 3.1070380314717982, -2.1861934003966104], [-2.2659422418234954, -1.026654842279342, -2.2798191986199234], [-4.202897212097876, -1.2125119740284835, -1.3132989952892569], [-2.6212054717929294, -2.6819879407832747, 0.09005605494850825], [-2.713377980741128, -0.7960712600074706, 1.1536489463716084], [-0.37104516182828085, -1.9490799840121111, 1.2553873021032873], [-0.24552532002049343, -2.1655580448419767, -0.8973245689258768], [1.577042443617164, 0.0, 0.0], [2.2927181468939164, 1.3915527243580552, 0.0], [3.782335574419713, 1.3186147352454594, -0.4807937764724515], [4.419592497958774, 0.34925073248439575, 0.14976318527134314], [3.8473329977166304, 1.1104432852892014, -1.7795740852228874], [4.368010411051578, 2.468434085335207, -0.20560554420872484], [2.2938967322202246, 1.8704189044736, 1.240568989312629], [1.6494649440008786, 2.235214894314341, -0.8029305726284841], [1.9974224573334871, -0.6906780683055302, 1.053572223549289], [1.9277183224308905, -0.6529932317206227, -1.1102241252095384], [-0.3501493572534786, 0.5705349971623116, -1.1530217920585761], [-0.42668432219275393, 0.7576153073313032, 1.0049834283127301], [-2.8416156026699833, 3.6456561923865074, -2.209782590839316]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0450', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
