import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0196'
logfile = 'conf/5009017845242299296281_0196.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863835, 0.7718203945763832, 1.163533622908847], [-0.3466020415139052, 0.2926212946843425, 2.625937724192379], [-0.9873929842445267, 1.0841115196807767, 3.819465912284932], [-2.462172495902625, 1.515549635617316, 3.5056638996836136], [-3.0963703180309756, 1.7977976497721433, 4.630222090259054], [-2.4735048074690122, 2.582029157264965, 2.7110614524348216], [-3.414321394728243, 0.15733163429363328, 2.6279853329745646], [-4.8025540278360275, 0.44673458194602744, 2.7351632135037294], [-2.747125159351999, -0.08691654408773691, 1.4030164425511547], [-3.0690607739072786, -1.0296540093275808, 3.6038966073004945], [-0.9969507321948253, 0.3041218762439119, 4.896738560267487], [-0.2788854647228145, 2.1786256533403403, 4.065404334434962], [-0.7784534903451429, -0.9699349115591107, 2.691387238785714], [0.9731582247379417, 0.30374941440734815, 2.768795082645789], [-2.0119541879597196, 0.605945531805917, 1.047594832227991], [-0.4099470658637787, 2.065478369993785, 1.0587059160250714], [1.5770424436171657, 0.0, 0.0], [2.292718146893917, 1.3915527243580545, 0.0], [2.341079856722327, 2.0598526928949283, 1.4165023767064702], [1.1453916410070388, 2.033172555891596, 1.9759511228307913], [3.202217992284645, 1.4422259961598272, 2.1984030628363262], [2.7215555591495084, 3.315205977928796, 1.2745358845394015], [1.6292323391939778, 2.2122558673105734, -0.8090479336198899], [3.545586830094382, 1.2600392214310687, -0.4280914688619813], [1.9974224573334838, -0.6906780683055262, 1.0535722235493015], [1.9277183224308974, -0.6529932317206231, -1.1102241252095306], [-0.3501493572534738, -1.2838136616209443, 0.08241309473865076], [-0.4266843221927563, 0.49153352563554553, -1.158605816601229], [-2.969506639336355, -1.8543408990487456, 3.1079703074253353]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0196', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
