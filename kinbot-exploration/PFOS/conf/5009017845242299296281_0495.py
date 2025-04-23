import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0495'
logfile = 'conf/5009017845242299296281_0495.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863867, -1.3935598728845953, 0.08664925740765313], [-2.270962283629197, -1.3932545648232102, 0.0536763686732784], [-2.969991788512713, -0.7382819684025391, -1.188880669680403], [-2.1956407956954163, -1.041846620920216, -2.518472411782727], [-1.1490862130831823, -0.2414432051524345, -2.6213306672491723], [-1.7886819047330238, -2.30789580125618, -2.5385715309252825], [-3.27647684277683, -0.7805268783930566, -4.030149878245057], [-4.15040995962661, -1.8959837637135948, -4.1502658724932155], [-3.7101390815568434, 0.5671491309313871, -3.997543004802992], [-2.1573635071723256, -0.8947422176778465, -5.132457833470466], [-4.2015349880113835, -1.2287045207786609, -1.2974100223711258], [-3.024435797927934, 0.5781431133077395, -1.0323857637828722], [-2.621205471792933, -2.681987940783274, 0.09005605494850584], [-2.7133779807411313, -0.7960712600074644, 1.1536489463716038], [-0.371045161828289, -1.94907998401211, 1.2553873021032915], [-0.24552532002049376, -2.1655580448419802, -0.897324568925876], [1.5770424436171646, 0.0, 0.0], [2.2927181468939195, 1.391552724358053, 0.0], [3.782335574419718, 1.3186147352454536, -0.48079377647243865], [4.419592497958777, 0.3492507324843841, 0.14976318527135118], [3.847332997716636, 1.1104432852891906, -1.7795740852228839], [4.368010411051582, 2.4684340853352014, -0.20560554420871946], [2.2938967322202197, 1.8704189044736035, 1.240568989312633], [1.6494649440008842, 2.235214894314334, -0.8029305726284856], [1.9974224573334798, -0.6906780683055287, 1.0535722235492981], [1.927718322430894, -0.652993231720626, -1.1102241252095308], [-0.3501493572534747, 0.5705349971623105, -1.1530217920585804], [-0.42668432219275515, 0.7576153073313063, 1.0049834283127268], [-2.1216872025558517, -1.795731004944968, -5.483041337366235]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0495', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
