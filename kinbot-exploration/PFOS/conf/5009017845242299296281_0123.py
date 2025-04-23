import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0123'
logfile = 'conf/5009017845242299296281_0123.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586388, -1.393559872884595, 0.08664925740765311], [-2.270962283629198, -1.3932545648232095, 0.05367636867327837], [-2.969991788512715, -0.738281968402538, -1.1888806696804026], [-2.1956407956954176, -1.041846620920214, -2.518472411782727], [-1.1490862130831847, -0.24144320515243145, -2.6213306672491714], [-1.788681904733024, -2.30789580125618, -2.5385715309252825], [-3.2764768427768316, -0.7805268783930518, -4.030149878245055], [-4.016406191958046, 0.4198604867198433, -3.845743425788588], [-2.4707732998027305, -1.0496883856983303, -5.163177086864583], [-4.245179114609112, -2.0064309472087296, -3.8318628327349145], [-4.201534988011386, -1.2287045207786584, -1.2974100223711231], [-3.0244357979279353, 0.5781431133077406, -1.0323857637828722], [-2.621205471792935, -2.6819879407832716, 0.09005605494850577], [-2.7133779807411322, -0.7960712600074628, 1.153648946371604], [-0.37104516182829006, -1.9490799840121107, 1.2553873021032897], [-0.24552532002049537, -2.1655580448419807, -0.8973245689258771], [1.5770424436171635, 0.0, 0.0], [2.2927181468939195, 1.3915527243580517, 0.0], [3.782335574419717, 1.3186147352454536, -0.48079377647243804], [4.419592497958777, 0.34925073248438177, 0.1497631852713489], [3.847332997716635, 1.1104432852891921, -1.7795740852228854], [4.36801041105158, 2.468434085335203, -0.20560554420871752], [2.2938967322202193, 1.8704189044736024, 1.2405689893126342], [1.6494649440008828, 2.235214894314335, -0.8029305726284834], [1.9974224573334805, -0.6906780683055295, 1.053572223549298], [1.9277183224308918, -0.6529932317206254, -1.110224125209532], [-0.3501493572534755, 0.5705349971623114, -1.1530217920585797], [-0.4266843221927548, 0.7576153073313063, 1.0049834283127268], [-5.051675770388328, -1.727201940695784, -3.3762789295570985]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0123', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
