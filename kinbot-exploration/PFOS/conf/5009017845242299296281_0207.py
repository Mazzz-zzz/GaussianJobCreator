import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0207'
logfile = 'conf/5009017845242299296281_0207.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863822, 0.6217394783082162, -1.2501828803164983], [-0.3976197158559575, -0.07566485901595582, -2.633910119820666], [1.0879414097563231, -0.0836723585219182, -3.138512306045849], [1.473042089945533, 1.2684473752977068, -3.833533780919367], [1.0004737759019138, 1.2895411745485015, -5.067525751424902], [0.9773040172258539, 2.2971743117836962, -3.1517239144922664], [3.332923574810669, 1.5027949558458802, -3.9206280107056743], [3.9198025913470134, 0.2649021047134447, -4.302236863950117], [3.5618854995970763, 2.7292422699378944, -4.590523670045682], [3.604488890927199, 1.7390700069393499, -2.387341015677571], [1.8947536492841386, -0.26661420225360427, -2.0970816341016425], [1.2598103705638535, -1.0641590715978342, -4.015744696739176], [-0.7843720924895164, -1.3445493546988763, -2.475221044816152], [-1.1581983763610557, 0.4930222734242936, -3.5614953002063556], [-0.2560445759534297, 1.8735740976390431, -1.3916701657561907], [-2.0076024771874463, 0.6427130616946742, -1.0543092166280577], [1.5770424436171657, 0.0, 0.0], [2.292718146893915, 1.3915527243580572, 0.0], [3.782335574419713, 1.3186147352454638, -0.48079377647244387], [4.419592497958774, 0.34925073248439853, 0.14976318527134896], [3.8473329977166326, 1.1104432852892019, -1.7795740852228776], [4.368010411051575, 2.4684340853352102, -0.2056055442087164], [2.2938967322202166, 1.870418904473609, 1.240568989312632], [1.6494649440008784, 2.23521489431434, -0.8029305726284856], [1.997422457333487, -0.690678068305523, 1.0535722235492986], [1.9277183224308965, -0.6529932317206211, -1.1102241252095326], [-0.35014935725347246, 0.7132786644586352, 1.0706086973199338], [-0.4266843221927524, -1.2491488329668536, 0.15362238828850272], [3.866405686235689, 0.9110912568991308, -1.9609454616280861]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0207', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
