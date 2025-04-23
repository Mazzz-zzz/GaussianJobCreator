import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0022'
logfile = 'conf/5009017845242299296281_0022.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863854, -1.3935598728845977, 0.0866492574076508], [-0.39761971585595896, -2.2432006455416373, 1.3824827499919154], [-1.1233882121466365, -3.6270084497455466, 1.5236046658601572], [-1.2235971086413773, -4.373754595344558, 0.14813618926815184], [-0.09347781129709948, -4.23332970793658, -0.5223958143256597], [-1.4669716767588572, -5.667130495747957, 0.3400935814215637], [-2.615554613648525, -3.6969298765090786, -0.9130075028596122], [-2.5918128170423977, -2.277701250039457, -0.8255006990495476], [-2.5907424079415757, -4.4171195142833835, -2.132051395317333], [-3.829529432267254, -4.234910948620465, -0.06602218089794128], [-0.4343595853052552, -4.389746980891084, 2.367678736664978], [-2.3533195674001606, -3.444774504497558, 1.9867347254034686], [0.9186702629217385, -2.470358014913694, 1.357831014730945], [-0.6863393184328452, -1.5018460017288915, 2.44518379143072], [-0.2560445759534338, -2.142008766053285, -0.9267276815498259], [-2.0076024771874494, -1.2344150958913078, -0.029451230457634058], [1.577042443617163, 0.0, 0.0], [2.292718146893914, 1.3915527243580543, 0.0], [1.6005215470082486, 2.440721904563898, -0.9357086002340341], [1.376069276137161, 1.9145350254105769, -2.1257143081021446], [0.4611709133627495, 2.8519467789832547, -0.4188289776134507], [2.404864907693506, 3.478080011182909, -1.0689303403306827], [3.540059223330461, 1.2295174124846042, -0.43152105569274757], [2.3088468039522527, 1.896094738758375, 1.231022041490469], [1.9974224573334824, -0.6906780683055298, 1.0535722235492924], [1.927718322430893, -0.6529932317206277, -1.1102241252095366], [-0.3501493572534765, 0.5705349971623109, -1.1530217920585841], [-0.42668432219275504, 0.7576153073313039, 1.0049834283127272], [-4.144389018860769, -3.5512782233811278, 0.5418200061934176]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0022', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
