import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0208'
logfile = 'conf/5009017845242299296281_0208.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863835, 0.7718203945763874, 1.1635336229088442], [-0.39761971585595535, 2.3188655045575954, 1.251427369828739], [1.0879414097563276, 2.759867566386749, 1.496793764948375], [1.4730420899455359, 2.685713952893136, 3.0152745408311827], [0.9614495666061569, 1.593462181345954, 3.5554236875418703], [2.795349477668079, 2.6726461892025335, 3.157515345382904], [0.8248630198382841, 4.1642430621827256, 3.9720811551337936], [0.8535738017776328, 3.844765015123364, 5.357555974918535], [1.4300365113903912, 5.314807706377585, 3.4105450331439475], [-0.6687885391477252, 4.120188233951884, 3.4744800877542907], [1.8947536492841413, 1.949433070068612, 0.8176461448894707], [1.2598103705638581, 4.009816458287683, 1.0862835586981867], [-0.7843720924895141, 2.815878982142089, 0.07319662459686575], [-1.1581983763610528, 2.837834268725452, 2.207717463520164], [-0.2560445759534298, 0.26843466841424896, 2.3183978473060125], [-2.0076024771874463, 0.5917020341966372, 1.0837604470856899], [1.5770424436171646, 0.0, 0.0], [2.2927181468939186, 1.3915527243580539, 0.0], [2.3410798567223323, 2.0598526928949292, 1.4165023767064728], [1.1453916410070413, 2.0331725558916, 1.9759511228307898], [3.20221799228465, 1.4422259961598276, 2.198403062836326], [2.7215555591495155, 3.315205977928793, 1.2745358845394017], [1.6292323391939854, 2.212255867310575, -0.8090479336198888], [3.545586830094386, 1.2600392214310652, -0.4280914688619793], [1.997422457333485, -0.690678068305526, 1.0535722235493006], [1.9277183224308971, -0.6529932317206302, -1.1102241252095295], [-0.3501493572534779, -1.2838136616209426, 0.08241309473865065], [-0.42668432219275754, 0.4915335256355425, -1.1586058166012303], [-1.2203053177951182, 3.6410898774678757, 4.108716947595097]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0208', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
