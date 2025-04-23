import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0415'
logfile = 'conf/5009017845242299296281_0415.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863822, 0.7718203945763782, 1.1635336229088529], [-0.346602041513899, 0.2926212946843335, 2.6259377241923825], [-0.9873929842445199, 1.084111519680764, 3.819465912284939], [-1.027435612376097, 0.2256586926736543, 5.1313502405463245], [-2.040447113782516, -0.6216591292156244, 5.082748287868311], [0.10871386352755073, -0.45017867605864714, 5.277323731228234], [-1.238879480957569, 1.2954728144902385, 6.658587146306561], [0.01404476657998734, 1.8990341837743374, 6.955783228465047], [-2.446767008943772, 2.0156461600603524, 6.49210938889213], [-1.507723043713296, 0.1487338720537606, 7.70417576831499], [-0.26149276674436317, 2.1741733763092608, 4.051469307654182], [-2.2308779150238935, 1.4310362312190683, 3.5131682254871475], [-0.7784534903451333, -0.9699349115591188, 2.6913872387857136], [0.9731582247379487, 0.30374941440734077, 2.76879508264579], [-2.011954187959718, 0.6059455318059114, 1.0475948322279984], [-0.40994706586378044, 2.065478369993781, 1.058705916025081], [1.5770424436171684, 0.0, 0.0], [2.292718146893911, 1.3915527243580617, 0.0], [2.3410798567223075, 2.0598526928949368, 1.4165023767064804], [1.1453916410070277, 2.0331725558915816, 1.9759511228307964], [3.202217992284636, 1.4422259961598334, 2.1984030628363307], [2.7215555591494756, 3.315205977928803, 1.2745358845394135], [1.629232339193957, 2.21225586731058, -0.8090479336198825], [3.545586830094372, 1.260039221431096, -0.42809146886198013], [1.997422457333489, -0.690678068305526, 1.0535722235492968], [1.9277183224308994, -0.6529932317206175, -1.110224125209532], [-0.3501493572534725, -1.2838136616209455, 0.08241309473864855], [-0.4266843221927586, 0.491533525635548, -1.1586058166012245], [-1.8591196186057128, -0.6355414723710527, 7.259884177153038]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0415', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
