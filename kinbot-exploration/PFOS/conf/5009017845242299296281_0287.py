import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0287'
logfile = 'conf/5009017845242299296281_0287.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863865, -1.3935598728845968, 0.08664925740765074], [-0.39761971585595896, -2.2432006455416373, 1.3824827499919154], [-1.1233882121466392, -3.627008449745546, 1.5236046658601567], [-1.223597108641382, -4.373754595344556, 0.14813618926815175], [-0.09347781129710513, -4.233329707936578, -0.5223958143256595], [-1.4669716767588634, -5.667130495747955, 0.3400935814215636], [-2.615554613648531, -3.696929876509074, -0.9130075028596076], [-2.399954397565101, -4.1213393355821575, -2.2530911441611727], [-3.818868641292418, -3.909618453221887, -0.19735696533495523], [-2.2605739139485053, -2.1663091364244758, -0.8046955336578572], [-0.4343595853052563, -4.389746980891084, 2.367678736664978], [-2.3533195674001606, -3.444774504497558, 1.9867347254034686], [0.9186702629217364, -2.4703580149136943, 1.3578310147309451], [-0.6863393184328432, -1.501846001728892, 2.4451837914307206], [-0.2560445759534373, -2.1420087660532836, -0.9267276815498264], [-2.0076024771874503, -1.2344150958913067, -0.029451230457632385], [1.577042443617162, 0.0, 0.0], [2.2927181468939137, 1.3915527243580539, 0.0], [2.3410798567223194, 2.0598526928949275, 1.4165023767064784], [1.1453916410070388, 2.0331725558915927, 1.9759511228307987], [3.2022179922846474, 1.4422259961598232, 2.198403062836332], [2.721555559149505, 3.315205977928795, 1.2745358845394041], [1.6292323391939731, 2.2122558673105743, -0.8090479336198788], [3.5455868300943765, 1.2600392214310765, -0.4280914688619866], [1.997422457333481, -0.6906780683055302, 1.0535722235492921], [1.9277183224308905, -0.6529932317206277, -1.110224125209538], [-0.3501493572534792, 0.5705349971623138, -1.1530217920585795], [-0.4266843221927526, 0.7576153073313049, 1.0049834283127286], [-1.7302689259614916, -1.8896913652933158, -1.5651043809936345]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0287', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
