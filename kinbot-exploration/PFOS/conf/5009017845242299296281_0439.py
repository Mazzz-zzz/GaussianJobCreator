import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0439'
logfile = 'conf/5009017845242299296281_0439.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863833, 0.7718203945763901, 1.163533622908846], [-2.2709622836291925, 0.7431123812655732, 1.1797556627388974], [-2.969991788512714, -0.6604598778102035, 1.2338112746327874], [-2.1956407956954216, -1.6601377768739842, 2.1615018464552516], [-2.991273105863037, -2.6553452534664084, 2.5122628029895897], [-1.1300133907845076, -2.146015754603099, 1.5312155105620204], [-1.584021457030936, -0.8268720901825432, 3.7277883499625974], [-0.43246898602355366, -0.05851116725545418, 3.402243565876866], [-2.723451699091648, -0.3073734435172929, 4.388981703860276], [-1.126754123236767, -2.11643388301138, 4.507839372479614], [-4.201534988011384, -0.5092377781085852, 1.7127943399246752], [-3.0244357979279393, -1.1831438545952289, 0.01550625874391525], [-2.6212054717929285, 1.4189848017416633, 2.2776416618875786], [-2.7133779807411287, 1.3971249246107003, 0.1125934612033508], [-0.3710451618282809, 2.0617372872159176, 1.0602591291106087], [-0.24552532002049118, 0.30567315029126746, 2.3240905646658536], [1.5770424436171644, 0.0, 0.0], [2.2927181468939195, 1.3915527243580503, 0.0], [3.782335574419718, 1.3186147352454527, -0.48079377647244614], [4.419592497958774, 0.3492507324843912, 0.1497631852713508], [3.847332997716642, 1.1104432852891868, -1.7795740852228759], [4.368010411051586, 2.468434085335196, -0.2056055442087213], [2.2938967322202175, 1.8704189044736053, 1.2405689893126297], [1.6494649440008855, 2.2352148943143337, -0.8029305726284919], [1.9974224573334802, -0.6906780683055262, 1.0535722235492981], [1.9277183224308905, -0.652993231720628, -1.1102241252095302], [-0.3501493572534797, -1.2838136616209412, 0.08241309473865056], [-0.4266843221927585, 0.4915335256355488, -1.1586058166012294], [-1.2998185286179962, -2.015702837755711, 5.454340484368947]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0439', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
