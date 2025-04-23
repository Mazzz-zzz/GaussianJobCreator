import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0133'
logfile = 'conf/5009017845242299296281_0133.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863826, -1.393559872884598, 0.0866492574076533], [-2.2709622836291934, -1.3932545648232142, 0.05367636867327484], [-3.0203184489305848, -0.6573324959322145, 1.2192077454993062], [-4.480058160568046, -0.2559865728509693, 0.8097366435926585], [-5.037241623014005, -1.2341925500908644, 0.11752164928389947], [-5.212262690872639, -0.006028967204993568, 1.891539657886009], [-4.5062276883420225, 1.2924765630661825, -0.25009372190774], [-3.4550410105005875, 1.1969601661214013, -1.2031489113135623], [-5.860046646976047, 1.5391792464343612, -0.5842625189400881], [-4.091173519573295, 2.332943248660112, 0.8570204289212924], [-3.086116569243988, -1.4709098648754981, 2.2693287026335534], [-2.363928874392459, 0.4475154150817301, 1.5491980808953385], [-2.596140664371263, -0.7866250332999407, -1.0913507571730163], [-2.6897770261666514, -2.652128594710337, 0.007528686185012738], [-0.3710451618282821, -1.9490799840121094, 1.2553873021032862], [-0.24552532002049005, -2.1655580448419807, -0.897324568925875], [1.577042443617166, 0.0, 0.0], [2.2927181468939173, 1.3915527243580554, 0.0], [3.782335574419715, 1.3186147352454636, -0.4807937764724377], [4.419592497958777, 0.34925073248439287, 0.1497631852713529], [3.8473329977166446, 1.110443285289192, -1.7795740852228716], [4.3680104110515785, 2.4684340853352063, -0.2056055442087097], [2.293896732220217, 1.8704189044736048, 1.240568989312635], [1.6494649440008846, 2.235214894314336, -0.8029305726284817], [1.9974224573334833, -0.6906780683055254, 1.0535722235493008], [1.927718322430896, -0.6529932317206258, -1.1102241252095315], [-0.3501493572534737, 0.5705349971623073, -1.1530217920585817], [-0.42668432219275565, 0.7576153073313021, 1.004983428312725], [-4.319940010896681, 1.9982282447497632, 1.7354247465910544]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0133', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
