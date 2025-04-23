import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0326'
logfile = 'conf/5009017845242299296281_0326.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586382, -1.3935598728845986, 0.08664925740765334], [-2.270962283629193, -1.3932545648232162, 0.05367636867328111], [-2.969991788512713, -0.738281968402546, -1.1888806696803988], [-3.0567398040787164, 0.821610279220505, -1.0517817396959952], [-4.057927071959132, 1.147154978605416, -0.25294092550737673], [-1.9208515428855046, 1.306979118885893, -0.5587553832825164], [-3.3438782542772705, 1.6524617020879715, -2.709762666803327], [-3.8032649569633916, 2.9768077460426903, -2.4699660424178815], [-2.232158454102361, 1.3342071169109517, -3.5271670199608556], [-4.559261845011894, 0.7765380181382833, -3.1957081957240323], [-2.2659422418234927, -1.0266548422793464, -2.279819198619925], [-4.2028972120978745, -1.21251197402849, -1.313298995289259], [-2.6212054717929245, -2.681987940783279, 0.09005605494850602], [-2.7133779807411265, -0.7960712600074752, 1.1536489463716086], [-0.3710451618282794, -1.9490799840121116, 1.2553873021032875], [-0.24552532002049038, -2.165558044841978, -0.8973245689258761], [1.577042443617164, 0.0, 0.0], [2.292718146893914, 1.3915527243580559, 0.0], [3.7823355744197102, 1.3186147352454671, -0.4807937764724479], [4.419592497958774, 0.3492507324843993, 0.1497631852713448], [3.84733299771663, 1.1104432852892032, -1.779574085222886], [4.368010411051574, 2.468434085335211, -0.20560554420872423], [2.2938967322202206, 1.8704189044736035, 1.24056898931263], [1.6494649440008753, 2.2352148943143426, -0.8029305726284834], [1.9974224573334873, -0.6906780683055258, 1.053572223549291], [1.9277183224308911, -0.6529932317206224, -1.1102241252095386], [-0.35014935725347857, 0.5705349971623119, -1.1530217920585757], [-0.42668432219275526, 0.7576153073313024, 1.004983428312729], [-4.54621353537201, -0.08733714268955665, -2.7603786253796616]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0326', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
