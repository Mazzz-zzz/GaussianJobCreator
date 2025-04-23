import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0214'
logfile = 'conf/5009017845242299296281_0214.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863817, 0.621739478308216, -1.250182880316502], [-2.2709622836291925, 0.6501421835576533, -1.2334320314121754], [-2.969991788512709, 1.3987418462127623, -0.04493060495238367], [-2.1956407956954127, 2.701984397794217, 0.35697056532747673], [-1.7814667182728683, 3.331102206613144, -0.728926556546902], [-2.9823571229144368, 3.508786314148332, 1.0634039512934506], [-0.6934923395346737, 2.322119367376902, 1.4157037259656826], [-1.1396204347114405, 2.042875586432666, 2.736848548196736], [0.13508964616168598, 1.4583142196598802, 0.6589677539086719], [-0.0343811651830556, 3.7524212609742937, 1.39695318395915], [-4.20153498801138, 1.7379422988872706, -0.41538431755354954], [-3.0244357979279344, 0.605000741287506, 1.0168795050389607], [-2.6212054717929267, 1.2630031390416303, -2.367697716836088], [-2.7133779807411282, -0.6010536646032223, -1.2662424075749594], [-0.3710451618282811, -0.11265730320379998, -2.315646431213901], [-0.2455253200204881, 1.8598848945507205, -1.426765995739979], [1.577042443617165, 0.0, 0.0], [2.2927181468939177, 1.3915527243580574, 0.0], [3.7823355744197156, 1.3186147352454578, -0.4807937764724498], [4.419592497958778, 0.34925073248439287, 0.14976318527134924], [3.847332997716637, 1.110443285289199, -1.7795740852228783], [4.368010411051581, 2.4684340853352063, -0.20560554420871868], [2.2938967322202197, 1.8704189044736055, 1.2405689893126315], [1.6494649440008822, 2.235214894314337, -0.8029305726284877], [1.9974224573334842, -0.6906780683055275, 1.0535722235492986], [1.9277183224308954, -0.6529932317206245, -1.1102241252095344], [-0.350149357253475, 0.7132786644586347, 1.0706086973199331], [-0.4266843221927575, -1.2491488329668505, 0.15362238828850008], [0.9300647801495169, 3.6778316151302453, 1.381371465170074]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0214', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
