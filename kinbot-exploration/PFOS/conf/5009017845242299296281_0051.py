import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0051'
logfile = 'conf/5009017845242299296281_0051.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863854, -1.3935598728845977, 0.0866492574076508], [-0.39761971585595896, -2.2432006455416373, 1.3824827499919154], [-1.1233882121466365, -3.6270084497455466, 1.5236046658601572], [-1.2235971086413773, -4.373754595344558, 0.14813618926815184], [-0.09347781129709948, -4.23332970793658, -0.5223958143256597], [-1.4669716767588572, -5.667130495747957, 0.3400935814215637], [-2.615554613648525, -3.6969298765090786, -0.9130075028596122], [-2.591812817042398, -2.2777012500394576, -0.8255006990495477], [-2.5907424079415757, -4.4171195142833835, -2.132051395317333], [-3.829529432267254, -4.234910948620465, -0.06602218089794128], [-0.4343595853052552, -4.389746980891084, 2.367678736664978], [-2.3533195674001606, -3.444774504497558, 1.9867347254034686], [0.9186702629217385, -2.470358014913694, 1.357831014730945], [-0.6863393184328452, -1.5018460017288915, 2.44518379143072], [-0.2560445759534338, -2.142008766053285, -0.9267276815498259], [-2.0076024771874494, -1.2344150958913078, -0.029451230457634058], [1.577042443617163, 0.0, 0.0], [2.292718146893914, 1.3915527243580543, 0.0], [3.7823355744197174, 1.3186147352454554, -0.48079377647244337], [4.419592497958774, 0.34925073248438976, 0.1497631852713488], [3.8473329977166393, 1.1104432852891948, -1.7795740852228805], [4.368010411051575, 2.4684340853352085, -0.20560554420871358], [2.293896732220216, 1.8704189044736026, 1.2405689893126322], [1.649464944000882, 2.2352148943143395, -0.802930572628486], [1.997422457333482, -0.6906780683055297, 1.0535722235492926], [1.9277183224308931, -0.6529932317206281, -1.1102241252095366], [-0.3501493572534765, 0.5705349971623109, -1.1530217920585841], [-0.42668432219275504, 0.7576153073313039, 1.0049834283127272], [-4.5588438490409295, -4.489804034275591, -0.6483385751951364]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0051', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
