import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0017'
logfile = 'conf/5009017845242299296281_0017.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863821, 0.771820394576383, 1.1635336229088513], [-0.3976197158559556, 2.318865504557591, 1.2514273698287515], [-0.7491833517666582, 3.196572017163505, -0.0006164121359253572], [0.37513150925003946, 3.1292034692660806, -1.0919288014421236], [-0.10815375808507582, 3.5174536458638235, -2.2590616532741867], [1.3965384052652352, 3.910843918740282, -0.7533048196106789], [1.0485429431935351, 1.3891449930549489, -1.2928687703755177], [1.7447219915669392, 1.3218845923418672, -2.53111017612474], [1.6092624577838988, 1.0187675041004514, -0.0463652953377202], [-0.3233323711935339, 0.6302082600302853, -1.4429213839823267], [-1.8820373528306364, 2.749164657044343, -0.5349009981858434], [-0.8963077255402097, 4.463477707111217, 0.3650804471736083], [-1.1305489784219276, 2.7568196252985233, 2.278967330580903], [0.8878652548597813, 2.4842258206717505, 1.5386269105814951], [-0.25604457595342656, 0.2684346684142412, 2.3183978473060174], [-2.007602477187446, 0.5917020341966347, 1.0837604470856967], [1.5770424436171644, 0.0, 0.0], [2.2927181468939155, 1.3915527243580548, 0.0], [3.7823355744197125, 1.3186147352454656, -0.48079377647244526], [4.419592497958776, 0.3492507324843961, 0.14976318527134902], [3.8473329977166335, 1.1104432852892008, -1.779574085222876], [4.368010411051578, 2.468434085335209, -0.20560554420871702], [2.293896732220216, 1.8704189044736053, 1.240568989312635], [1.6494649440008784, 2.2352148943143404, -0.8029305726284862], [1.9974224573334824, -0.690678068305528, 1.053572223549298], [1.927718322430895, -0.6529932317206233, -1.1102241252095328], [-0.3501493572534766, -1.2838136616209455, 0.08241309473865198], [-0.42668432219276037, 0.49153352563554625, -1.1586058166012247], [-0.2744489890056814, -0.23898830032312024, -1.0209378047855933]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0017', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
