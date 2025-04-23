import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0407'
logfile = 'conf/5009017845242299296281_0407.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863851, 0.6217394783082126, -1.2501828803165012], [-0.3976197158559567, -0.07566485901595647, -2.6339101198206682], [-1.1233882121466345, 0.49402387891336585, -3.9028837901505384], [-0.3678111583333635, 0.11496132015469422, -5.223868290930238], [0.6688104135818967, 0.9167943976383535, -5.394107462597693], [0.05149838568730978, -1.1461400629926168, -5.171879366435433], [-1.475740760567105, 0.27151298615619196, -6.730403024509696], [-2.224617795791914, 1.4743444314632588, -6.608703771827486], [-0.6859569336734569, -0.06064640797263146, -7.857908419315947], [-2.428681401794962, -0.9483036760987582, -6.4397243796276955], [-2.351535285418466, -0.013330338437770943, -3.9597873669787513], [-1.1885301103712396, 1.8171105205812563, -3.827863304783673], [0.9186702629217409, 0.05926285465343763, -2.8183083047232307], [-0.686339318432839, -1.3666682794365137, -2.5232286857846646], [-0.25604457595343316, 1.8735740976390431, -1.3916701657561947], [-2.007602477187449, 0.642713061694672, -1.0543092166280623], [1.5770424436171657, 0.0, 0.0], [2.29271814689391, 1.3915527243580619, 0.0], [1.6005215470082381, 2.440721904563901, -0.93570860023403], [1.3760692761371511, 1.9145350254105629, -2.1257143081021503], [0.4611709133627373, 2.8519467789832462, -0.41882897761345295], [2.404864907693479, 3.478080011182914, -1.0689303403306911], [3.5400592233304593, 1.2295174124846202, -0.43152105569275073], [2.308846803952244, 1.8960947387583846, 1.2310220414904678], [1.9974224573334896, -0.6906780683055197, 1.0535722235492988], [1.9277183224308956, -0.6529932317206208, -1.110224125209531], [-0.3501493572534747, 0.7132786644586332, 1.0706086973199331], [-0.42668432219275126, -1.2491488329668525, 0.15362238828850258], [-2.6831538097969703, -1.3805454738806904, -7.266992902244088]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0407', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
