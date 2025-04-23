import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0400'
logfile = 'conf/5009017845242299296281_0400.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863804, 0.7718203945763885, 1.1635336229088438], [-0.34660204151390467, 0.2926212946843473, 2.625937724192379], [1.1624535490467056, 0.3084161010830978, 3.054984450098302], [1.9172819807099757, 1.5569601130304869, 2.479695785178603], [2.2272170921886336, 1.3494186210778614, 1.211875333037296], [1.1566135117755312, 2.643560915797123, 2.5775845810075273], [3.5062602700141086, 1.9036546628801005, 3.4159725692225535], [3.169361687962342, 2.5315200986645476, 4.64671029423496], [4.306971521358234, 0.7402973694523728, 3.311962610817348], [4.110575342315997, 2.9867518250503053, 2.445278765545708], [1.2315520144708494, 0.3493659548263259, 4.38259134374583], [1.7685418932081407, -0.7841917084166115, 2.608559881799362], [-1.0117971004258965, 1.1239736000383767, 3.4329051118380067], [-0.8215282967123825, -0.9363359387499598, 2.788444311143955], [-2.0119541879597156, 0.6059455318059208, 1.0475948322279889], [-0.4099470658637755, 2.0654783699937878, 1.0587059160250727], [1.5770424436171657, 0.0, 0.0], [2.2927181468939173, 1.3915527243580552, 0.0], [3.782335574419716, 1.3186147352454594, -0.4807937764724469], [4.419592497958776, 0.3492507324843901, 0.14976318527135302], [3.8473329977166406, 1.110443285289195, -1.779574085222876], [4.36801041105158, 2.468434085335204, -0.20560554420871746], [2.293896732220217, 1.870418904473607, 1.240568989312631], [1.6494649440008877, 2.235214894314335, -0.8029305726284892], [1.9974224573334816, -0.6906780683055262, 1.0535722235493006], [1.927718322430899, -0.6529932317206311, -1.11022412520953], [-0.35014935725347374, -1.2838136616209455, 0.08241309473865083], [-0.42668432219275476, 0.49153352563554453, -1.158605816601228], [5.073276023416751, 2.9006779780220864, 2.4033135637482363]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0400', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
