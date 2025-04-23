import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0159'
logfile = 'conf/5009017845242299296281_0159.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863833, 0.7718203945763831, 1.163533622908849], [-0.3466020415139029, 0.29262129468434256, 2.6259377241923794], [1.1624535490467058, 0.30841610108309514, 3.0549844500983023], [1.3229555267971613, 0.3508819546540548, 4.6144782543247524], [2.5379668391891173, -0.044715161615564956, 4.951750214967398], [1.1168882370515762, 1.5841606881125143, 5.067755077302737], [0.08826989713528051, -0.7713860003073009, 5.473389986860927], [0.5381926157038321, -0.9914411247831094, 6.804395733399166], [-1.2025433687641451, -0.2930625394278117, 5.141178381105891], [0.34402262312628656, -2.0775804215303504, 4.631374959131036], [1.750051457246916, -0.7950580504033852, 2.6010405683398576], [1.7633282149005842, 1.376908489527518, 2.5473394429719383], [-1.0117971004258974, 1.1239736000383704, 3.432905111838009], [-0.821528296712381, -0.9363359387499647, 2.7884443111439534], [-2.011954187959717, 0.6059455318059157, 1.047594832227993], [-0.4099470658637775, 2.0654783699937855, 1.0587059160250751], [1.577042443617167, 0.0, 0.0], [2.292718146893914, 1.3915527243580588, 0.0], [3.7823355744197125, 1.3186147352454685, -0.4807937764724453], [4.419592497958776, 0.349250732484398, 0.14976318527135218], [3.847332997716639, 1.1104432852892019, -1.7795740852228785], [4.368010411051578, 2.4684340853352094, -0.20560554420871824], [2.2938967322202135, 1.8704189044736095, 1.2405689893126333], [1.649464944000882, 2.235214894314339, -0.8029305726284865], [1.9974224573334851, -0.6906780683055245, 1.0535722235492988], [1.9277183224308962, -0.6529932317206281, -1.1102241252095317], [-0.35014935725347296, -1.2838136616209455, 0.08241309473865083], [-0.4266843221927547, 0.49153352563554603, -1.1586058166012285], [0.7028940834630952, -1.8519322129328848, 3.7617448379826945]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0159', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
