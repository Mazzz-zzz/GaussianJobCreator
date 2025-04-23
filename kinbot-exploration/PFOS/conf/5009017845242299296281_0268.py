import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0268'
logfile = 'conf/5009017845242299296281_0268.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863814, -1.393559872884599, 0.0866492574076484], [-2.270962283629191, -1.3932545648232157, 0.05367636867327118], [-3.0203184489305848, -0.6573324959322174, 1.219207745499301], [-2.2823447039971083, 0.6594581441991485, 1.6446984279510883], [-1.8454205729265496, 1.2959521938724377, 0.5720966511971368], [-3.105173856347475, 1.4520556038051027, 2.325610678994922], [-0.8077838749090355, 0.3063105249128092, 2.7504201160021453], [-0.10865718900749041, -0.810908935992189, 2.2161228569361158], [-0.20471546093758505, 1.5538604817781554, 3.042454985227118], [-1.592498989935659, -0.15632231755508136, 4.035229701554743], [-4.245418844080691, -0.34008637241846496, 0.8099700744271577], [-3.094058093800262, -1.4520145440421857, 2.2791477007209346], [-2.596140664371263, -0.7866250332999424, -1.0913507571730208], [-2.6897770261666483, -2.6521285947103386, 0.007528686185005677], [-0.37104516182828196, -1.9490799840121096, 1.255387302103284], [-0.24552532002048594, -2.1655580448419802, -0.8973245689258781], [1.5770424436171657, 0.0, 0.0], [2.2927181468939146, 1.3915527243580608, 0.0], [3.7823355744197142, 1.3186147352454654, -0.4807937764724397], [4.419592497958777, 0.3492507324843972, 0.1497631852713539], [3.847332997716641, 1.110443285289201, -1.7795740852228725], [4.368010411051576, 2.468434085335217, -0.2056055442087087], [2.293896732220216, 1.8704189044736066, 1.2405689893126353], [1.6494649440008802, 2.2352148943143444, -0.8029305726284883], [1.9974224573334818, -0.6906780683055223, 1.0535722235493015], [1.9277183224308962, -0.6529932317206228, -1.1102241252095326], [-0.3501493572534757, 0.57053499716231, -1.1530217920585835], [-0.4266843221927585, 0.7576153073313048, 1.0049834283127228], [-1.1311056821209438, 0.13530928583754323, 4.833997258546415]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0268', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
