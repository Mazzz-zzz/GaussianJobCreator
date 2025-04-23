import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0342'
logfile = 'conf/5009017845242299296281_0342.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863828, 0.7718203945763873, 1.1635336229088467], [-2.270962283629193, 0.7431123812655708, 1.1797556627388979], [-2.969991788512713, -0.6604598778102062, 1.2338112746327872], [-2.1956407956954176, -1.6601377768739882, 2.1615018464552525], [-2.9912731058630335, -2.655345253466412, 2.512262802989589], [-1.1300133907845047, -2.1460157546031, 1.531215510562021], [-1.584021457030934, -0.8268720901825498, 3.7277883499625966], [-1.2617908272916074, -1.8390626909337808, 4.673340366944559], [-0.6920488892418222, 0.19745524149161767, 3.327169793796012], [-2.9335170727851736, -0.14092262343103906, 4.162388789369187], [-4.2015349880113835, -0.5092377781085887, 1.7127943399246761], [-3.0244357979279366, -1.1831438545952324, 0.015506258743914246], [-2.621205471792927, 1.4189848017416604, 2.2776416618875786], [-2.713377980741127, 1.3971249246106994, 0.11259346120335136], [-0.3710451618282802, 2.0617372872159154, 1.0602591291106098], [-0.24552532002049016, 0.3056731502912654, 2.324090564665854], [1.5770424436171655, 0.0, 0.0], [2.2927181468939186, 1.3915527243580526, 0.0], [2.3410798567223274, 2.059852692894926, 1.4165023767064748], [1.1453916410070426, 2.0331725558915967, 1.9759511228307907], [3.202217992284644, 1.4422259961598223, 2.198403062836328], [2.7215555591495137, 3.315205977928792, 1.2745358845394017], [1.6292323391939836, 2.2122558673105743, -0.8090479336198887], [3.5455868300943862, 1.2600392214310658, -0.4280914688619826], [1.9974224573334818, -0.6906780683055278, 1.0535722235492995], [1.927718322430894, -0.6529932317206304, -1.1102241252095295], [-0.3501493572534773, -1.2838136616209426, 0.08241309473865065], [-0.4266843221927564, 0.49153352563554675, -1.158605816601229], [-3.505752722677643, -0.008083531813221039, 3.393712964702961]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0342', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
