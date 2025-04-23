import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0357'
logfile = 'conf/5009017845242299296281_0357.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863848, -1.3935598728845966, 0.08664925740764574], [-2.2709622836291934, -1.3932545648232142, 0.0536763686732736], [-2.9970239643019028, -2.7837283183799224, 0.08524408716292907], [-2.2704853996782592, -3.7973876868622556, 1.0361497385314615], [-1.2136023415164636, -4.305894295055957, 0.42711675595095866], [-1.8869755643483324, -3.189249965050728, 2.1550727734485697], [-3.393595145601963, -5.2170258152991895, 1.5310390332371289], [-4.290607532258395, -4.7472249939132265, 2.5295814580844684], [-3.8019428663004198, -5.854656255566395, 0.3343618285185628], [-2.3119855550805224, -6.134911195652059, 2.215253730562105], [-4.239914773111564, -2.6098357638694933, 0.5256704367019929], [-3.023054582933855, -3.305441500441067, -1.1344147028460934], [-2.6427825570536165, -0.7109480038655843, 1.14043067756133], [-2.667789377892035, -0.7313203152515799, -1.02644886795114], [-0.3710451618282825, -1.949079984012116, 1.2553873021032806], [-0.24552532002049052, -2.165558044841976, -0.8973245689258843], [1.5770424436171642, 0.0, 0.0], [2.2927181468939173, 1.3915527243580539, 0.0], [1.600521547008264, 2.440721904563901, -0.9357086002340267], [1.3760692761371733, 1.9145350254105804, -2.1257143081021406], [0.4611709133627657, 2.8519467789832618, -0.4188289776134434], [2.4048649076935185, 3.4780800111829073, -1.068930340330693], [3.5400592233304677, 1.2295174124845962, -0.4315210556927503], [2.3088468039522625, 1.8960947387583689, 1.231022041490469], [1.9974224573334811, -0.6906780683055278, 1.053572223549295], [1.9277183224308923, -0.6529932317206169, -1.1102241252095337], [-0.35014935725347457, 0.5705349971623144, -1.1530217920585768], [-0.4266843221927555, 0.7576153073313006, 1.0049834283127321], [-2.5077751924305574, -7.067360936443369, 2.047425587807206]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0357', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
