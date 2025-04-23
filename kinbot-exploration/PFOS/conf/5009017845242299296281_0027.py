import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0027'
logfile = 'conf/5009017845242299296281_0027.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863848, 0.6217394783082156, -1.2501828803164987], [-0.346602041513907, 2.1278181305643273, -1.566386336981117], [-0.9873929842445301, 2.7656987490870772, -2.8486010727213635], [-0.1847769512407539, 2.396798184547748, -4.14448109650823], [0.1883244129367716, 1.1299271526239631, -4.0961581868459], [-0.9332368412513781, 2.597205265795233, -5.225473619786791], [1.3577185656456106, 3.448935518528457, -4.332571345974898], [2.2189737591343435, 2.8123260914604726, -5.268198909093286], [0.9276165596370575, 4.794030109958247, -4.437878402631883], [1.9474063794356815, 3.248787107693891, -2.8859367944322867], [-2.23050689855434, 2.3112475123060565, -2.9794957706412695], [-1.0036314587762374, 4.0874032189906675, -2.733434203914622], [-0.7784534903451426, 2.815777175989239, -0.505705345965246], [0.9731582247379391, 2.2459721722410166, -1.6474522505842988], [-2.0119541879597183, 0.6042709716797791, -1.0485616399675903], [-0.4099470658637772, -0.11587296658229848, -2.3181096973944286], [1.5770424436171644, 0.0, 0.0], [2.292718146893914, 1.3915527243580579, 0.0], [3.7823355744197125, 1.3186147352454631, -0.4807937764724484], [4.419592497958772, 0.34925073248439464, 0.14976318527134497], [3.8473329977166326, 1.1104432852892092, -1.7795740852228783], [4.368010411051575, 2.468434085335207, -0.20560554420871696], [2.2938967322202144, 1.8704189044736044, 1.2405689893126386], [1.649464944000878, 2.2352148943143426, -0.8029305726284831], [1.9974224573334873, -0.69067806830553, 1.0535722235492926], [1.9277183224308951, -0.6529932317206212, -1.1102241252095357], [-0.3501493572534749, 0.7132786644586303, 1.0706086973199371], [-0.42668432219275526, -1.2491488329668545, 0.1536223882884972], [1.2419200179155014, 3.0225452525993584, -2.2637905066442756]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0027', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
