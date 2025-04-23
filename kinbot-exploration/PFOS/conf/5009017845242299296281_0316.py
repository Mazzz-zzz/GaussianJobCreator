import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0316'
logfile = 'conf/5009017845242299296281_0316.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863844, -1.3935598728845988, 0.08664925740764715], [-0.39761971585595707, -2.2432006455416382, 1.3824827499919135], [-1.1233882121466345, -3.6270084497455497, 1.5236046658601508], [-0.3678111583333658, -4.58148330604694, 2.5123747217585426], [-1.1821704538825695, -5.537425913369774, 2.9240719555777956], [0.6879514164795761, -5.126087444421055, 1.9143180700580986], [0.26055307355107626, -3.6639458561435103, 4.023978754310711], [0.563863206493521, -4.621128198365776, 5.031071567710913], [1.1717646888843918, -2.683903805254325, 3.56076417960632], [-1.075302365678474, -2.9264726945328867, 4.414063631220663], [-2.351535285418468, -3.422611284169414, 1.9914380952175164], [-1.1885301103712396, -4.223582124447546, 0.340267780084486], [0.9186702629217401, -2.470358014913695, 1.3578310147309456], [-0.6863393184328435, -1.5018460017288928, 2.4451837914307175], [-0.25604457595343205, -2.1420087660532823, -0.9267276815498281], [-2.0076024771874477, -1.2344150958913078, -0.029451230457637898], [1.5770424436171655, 0.0, 0.0], [2.2927181468939173, 1.3915527243580557, 0.0], [3.782335574419715, 1.3186147352454598, -0.4807937764724466], [4.419592497958778, 0.3492507324843941, 0.14976318527134735], [3.8473329977166366, 1.1104432852891992, -1.7795740852228807], [4.368010411051581, 2.468434085335206, -0.2056055442087119], [2.2938967322202197, 1.870418904473599, 1.2405689893126286], [1.649464944000881, 2.2352148943143404, -0.8029305726284836], [1.9974224573334836, -0.6906780683055289, 1.053572223549294], [1.927718322430894, -0.6529932317206248, -1.1102241252095386], [-0.35014935725347845, 0.5705349971623136, -1.1530217920585817], [-0.42668432219275365, 0.7576153073313037, 1.0049834283127288], [-0.8827468948257032, -2.0460599667871278, 4.765856476747604]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0316', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
