import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0258'
logfile = 'conf/5009017845242299296281_0258.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863852, 0.6217394783082163, -1.2501828803165014], [-0.3466020415139046, 2.1278181305643256, -1.5663863369811188], [1.1624535490467054, 2.4914860914100188, -1.7945884035232567], [1.322955526797161, 3.820814416129081, -2.6111118136223204], [2.5379668391891146, 4.31069906016461, -2.437150641590285], [1.1168882370515751, 3.5967242930454923, -3.905800938233436], [0.08826989713527926, 5.125787773594593, -2.068655121040657], [-1.1777533685614712, 4.800388691205312, -2.628772551321797], [0.3103549562009413, 5.353582668209955, -0.6887477023436671], [0.68971223898328, 6.349584395209911, -2.8567741273869904], [1.7500514572469157, 2.650096233657924, -0.6119798150372631], [1.7633282149005813, 1.5176064249120464, -2.4661074521032575], [-1.0117971004258985, 2.4109962356139936, -2.689842246735281], [-0.8215282967123829, 2.8830315798638457, -0.5833314461381525], [-2.0119541879597196, 0.6042709716797823, -1.0485616399675937], [-0.40994706586377916, -0.11587296658229962, -2.318109697394431], [1.5770424436171668, 0.0, 0.0], [2.2927181468939155, 1.391552724358055, 0.0], [2.3410798567223234, 2.0598526928949283, 1.416502376706473], [1.145391641007036, 2.0331725558915887, 1.9759511228307924], [3.202217992284647, 1.4422259961598298, 2.198403062836323], [2.7215555591494986, 3.3152059779287946, 1.2745358845394068], [1.6292323391939714, 2.212255867310573, -0.8090479336198872], [3.54558683009438, 1.260039221431076, -0.42809146886198424], [1.9974224573334836, -0.6906780683055291, 1.0535722235492972], [1.9277183224308923, -0.6529932317206242, -1.1102241252095362], [-0.3501493572534752, 0.7132786644586331, 1.070608697319935], [-0.42668432219275654, -1.2491488329668528, 0.15362238828849809], [0.23162371590790237, 6.463356731815811, -3.701269922015864]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0258', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
