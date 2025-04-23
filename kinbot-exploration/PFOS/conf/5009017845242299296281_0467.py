import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0467'
logfile = 'conf/5009017845242299296281_0467.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863833, 0.7718203945763901, 1.163533622908846], [-2.2709622836291925, 0.7431123812655732, 1.1797556627388974], [-2.969991788512714, -0.6604598778102035, 1.2338112746327874], [-2.1956407956954216, -1.6601377768739842, 2.1615018464552516], [-2.991273105863037, -2.6553452534664084, 2.5122628029895897], [-1.1300133907845076, -2.146015754603099, 1.5312155105620204], [-1.584021457030936, -0.8268720901825432, 3.7277883499625974], [-0.4324689860235536, -0.058511167255454175, 3.4022435658768657], [-2.723451699091648, -0.3073734435172929, 4.388981703860276], [-1.126754123236767, -2.11643388301138, 4.507839372479614], [-4.201534988011384, -0.5092377781085852, 1.7127943399246752], [-3.0244357979279393, -1.1831438545952289, 0.01550625874391525], [-2.6212054717929285, 1.4189848017416633, 2.2776416618875786], [-2.7133779807411287, 1.3971249246107003, 0.1125934612033508], [-0.3710451618282809, 2.0617372872159176, 1.0602591291106087], [-0.24552532002049118, 0.30567315029126746, 2.3240905646658536], [1.5770424436171644, 0.0, 0.0], [2.2927181468939195, 1.3915527243580503, 0.0], [1.60052154700826, 2.4407219045638993, -0.9357086002340265], [1.3760692761371778, 1.9145350254105828, -2.125714308102138], [0.46117091336276483, 2.8519467789832555, -0.4188289776134505], [2.4048649076935122, 3.478080011182908, -1.0689303403306802], [3.540059223330469, 1.2295174124845976, -0.43152105569273935], [2.308846803952253, 1.896094738758368, 1.231022041490474], [1.9974224573334805, -0.6906780683055266, 1.0535722235492981], [1.9277183224308903, -0.652993231720628, -1.1102241252095304], [-0.3501493572534797, -1.2838136616209412, 0.08241309473865056], [-0.4266843221927585, 0.4915335256355488, -1.1586058166012294], [-0.17951094606808532, -2.2673221150183034, 4.381643130676575]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0467', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
