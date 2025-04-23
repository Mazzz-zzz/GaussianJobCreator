import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0004'
logfile = 'conf/5009017845242299296281_0004.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863841, 0.7718203945763799, 1.163533622908851], [-0.34660204151390284, 0.29262129468433784, 2.625937724192379], [1.1624535490467065, 0.3084161010830897, 3.054984450098303], [1.924360652536601, -0.9650386511798231, 2.5476426337726146], [1.6581688080331816, -1.9914892152919703, 3.336339714045335], [1.5600180827406167, -1.2571338783423756, 1.3023174109266629], [3.7832317877839663, -0.707808688603004, 2.5379379950928667], [4.116829961099136, 0.08349590716172445, 1.4043682780959856], [4.161464969515623, -0.39137245480168664, 3.8654034529571395], [4.217392380657199, -2.196099018049155, 2.2603788037803687], [1.7462981079693096, 1.383924747436085, 2.534025492845024], [1.2535978601361863, 0.34611839988489745, 4.378125052116879], [-1.011797100425898, 1.1239736000383642, 3.4329051118380103], [-0.8215282967123789, -0.9363359387499721, 2.7884443111439507], [-2.0119541879597187, 0.6059455318059094, 1.0475948322279949], [-0.40994706586378143, 2.06547836999378, 1.0587059160250756], [1.577042443617166, 0.0, 0.0], [2.292718146893911, 1.3915527243580585, 0.0], [3.7823355744197107, 1.31861473524547, -0.480793776472447], [4.4195924979587735, 0.34925073248439986, 0.1497631852713493], [3.8473329977166353, 1.1104432852892077, -1.779574085222878], [4.368010411051573, 2.468434085335214, -0.20560554420871752], [2.2938967322202113, 1.8704189044736086, 1.240568989312635], [1.649464944000876, 2.2352148943143395, -0.8029305726284833], [1.9974224573334853, -0.6906780683055245, 1.0535722235492966], [1.9277183224308958, -0.6529932317206225, -1.110224125209534], [-0.3501493572534724, -1.2838136616209472, 0.08241309473864752], [-0.4266843221927596, 0.49153352563554487, -1.1586058166012259], [3.527121769016021, -2.80763215765937, 2.5528101436859623]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0004', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
