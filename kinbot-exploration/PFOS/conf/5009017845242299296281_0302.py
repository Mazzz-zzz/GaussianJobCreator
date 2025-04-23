import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0302'
logfile = 'conf/5009017845242299296281_0302.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863848, -1.3935598728845966, 0.08664925740764574], [-2.2709622836291934, -1.3932545648232142, 0.0536763686732736], [-2.9970239643019028, -2.7837283183799224, 0.08524408716292907], [-2.2704853996782592, -3.7973876868622556, 1.0361497385314615], [-1.2136023415164636, -4.305894295055957, 0.42711675595095866], [-1.8869755643483324, -3.189249965050728, 2.1550727734485697], [-3.393595145601963, -5.2170258152991895, 1.5310390332371289], [-4.290607532258395, -4.7472249939132265, 2.5295814580844684], [-3.8019428663004198, -5.854656255566395, 0.3343618285185628], [-2.3119855550805224, -6.134911195652059, 2.215253730562105], [-4.239914773111564, -2.6098357638694933, 0.5256704367019929], [-3.023054582933855, -3.305441500441067, -1.1344147028460934], [-2.6427825570536165, -0.7109480038655843, 1.14043067756133], [-2.667789377892035, -0.7313203152515799, -1.02644886795114], [-0.3710451618282825, -1.949079984012116, 1.2553873021032806], [-0.24552532002049052, -2.165558044841976, -0.8973245689258843], [1.5770424436171642, 0.0, 0.0], [2.2927181468939173, 1.3915527243580539, 0.0], [3.782335574419714, 1.3186147352454642, -0.48079377647243954], [4.4195924979587735, 0.34925073248439475, 0.14976318527135252], [3.8473329977166366, 1.1104432852892048, -1.7795740852228772], [4.368010411051586, 2.4684340853352045, -0.20560554420870297], [2.2938967322202184, 1.870418904473599, 1.2405689893126386], [1.6494649440008815, 2.23521489431434, -0.8029305726284788], [1.9974224573334807, -0.6906780683055276, 1.0535722235492955], [1.9277183224308918, -0.6529932317206169, -1.110224125209534], [-0.35014935725347457, 0.5705349971623144, -1.1530217920585768], [-0.4266843221927555, 0.7576153073313006, 1.0049834283127321], [-1.4323916917877726, -5.936312815526436, 1.864774479014726]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0302', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
