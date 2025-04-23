import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0187'
logfile = 'conf/5009017845242299296281_0187.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863828, 0.6217394783082129, -1.2501828803165027], [-0.3466020415139064, 2.1278181305643233, -1.5663863369811226], [-0.7363023803695226, 3.202618946498001, -0.4919021834050389], [0.11232051632788868, 4.512501111822521, -0.6456326331137635], [0.2504583136006516, 4.810364670493019, -1.9257487186535769], [-0.4794659588750192, 5.524080936392357, -0.016828158858703337], [1.826573903151089, 4.323642832150675, 0.09418572775568974], [2.3338654001233086, 3.0484494962787845, -0.2785987958737094], [2.5077972431396933, 5.544793764029225, -0.12960420560801503], [1.4148838370527916, 4.278872665623413, 1.613741436334388], [-2.0209682759934697, 3.5132134413207536, -0.6407516003085895], [-0.5312999584009223, 2.716571178249477, 0.7255383184087602], [0.9809271994195052, 2.1566775839163785, -1.713562194923179], [-0.9145558447621731, 2.4557522019201086, -2.720590366825869], [-2.011954187959719, 0.6042709716797771, -1.048561639967592], [-0.4099470658637758, -0.11587296658230434, -2.318109697394432], [1.5770424436171633, 0.0, 0.0], [2.292718146893912, 1.3915527243580588, 0.0], [3.7823355744197116, 1.3186147352454651, -0.4807937764724501], [4.419592497958771, 0.3492507324843992, 0.14976318527135185], [3.847332997716634, 1.1104432852892046, -1.7795740852228805], [4.368010411051575, 2.4684340853352116, -0.20560554420871802], [2.2938967322202135, 1.8704189044736088, 1.2405689893126328], [1.6494649440008775, 2.2352148943143386, -0.802930572628489], [1.9974224573334833, -0.6906780683055235, 1.0535722235492966], [1.927718322430894, -0.652993231720623, -1.11022412520953], [-0.3501493572534777, 0.7132786644586345, 1.070608697319933], [-0.4266843221927575, -1.2491488329668519, 0.15362238828849686], [2.0818728570573652, 4.722956797568944, 2.1558459628267577]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0187', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
