import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0001'
logfile = 'conf/5009017845242299296281_0001.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863854, 0.7718203945763867, 1.1635336229088458], [-0.39761971585595896, 2.3188655045575954, 1.251427369828739], [-1.1233882121466365, 3.1329845708321837, 2.379279124290362], [-1.2235971086413773, 2.315167000798339, 3.713714494853226], [-0.09347781129709946, 1.6642568079316327, 3.9273689768312567], [-1.466971676758857, 3.1280949290491105, 4.737832185168432], [-2.615554613648525, 1.0577772469323425, 3.65813894049633], [-2.5918128170423977, 0.4239460488010215, 2.385297494290513], [-2.5907424079415757, 0.362149086622849, 4.891363408580055], [-3.829529432267254, 2.060278588439385, 3.700551554719137], [-0.4343595853052552, 4.2453434243976735, 2.6177930333052086], [-2.3533195674001606, 3.442949995028893, 1.9898948685020872], [0.9186702629217384, 2.411095160260257, 1.460477289992268], [-0.6863393184328452, 2.8685142811654, 0.07804489435393196], [-0.2560445759534338, 0.2684346684142513, 2.3183978473060147], [-2.0076024771874494, 0.5917020341966377, 1.0837604470856907], [1.577042443617163, 0.0, 0.0], [2.292718146893916, 1.3915527243580534, 0.0], [3.7823355744197125, 1.3186147352454558, -0.4807937764724462], [4.419592497958775, 0.3492507324843954, 0.14976318527135218], [3.8473329977166375, 1.1104432852891866, -1.7795740852228763], [4.368010411051574, 2.4684340853352085, -0.20560554420872446], [2.293896732220214, 1.8704189044736101, 1.2405689893126275], [1.649464944000882, 2.235214894314335, -0.8029305726284922], [1.9974224573334818, -0.6906780683055261, 1.0535722235493066], [1.9277183224308945, -0.6529932317206307, -1.110224125209528], [-0.35014935725347796, -1.283813661620946, 0.0824130947386566], [-0.42668432219275726, 0.49153352563553976, -1.1586058166012296], [-4.144389018860769, 2.2448690013327184, 2.8045871542578062]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0001', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
