import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0292'
logfile = 'conf/5009017845242299296281_0292.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863817, 0.6217394783082073, -1.250182880316503], [-0.3466020415139008, 2.1278181305643185, -1.5663863369811326], [-0.736302380369517, 3.202618946498001, -0.49190218340505193], [-0.5090926255515532, 2.6747258360153725, 0.9673064206754407], [-1.5087734085102669, 1.8856137906181645, 1.3200704424649146], [0.6350825859976149, 2.0009100837141847, 1.043039225853966], [-0.41185458969835453, 4.078131100526796, 2.2093457755430825], [-1.4318479285867929, 5.018976442093455, 1.8981309949083491], [-0.2399351960211971, 3.487262819072511, 3.4847575395224424], [0.977471611715914, 4.675425924366624, 1.7693977907915226], [0.018370741190138654, 4.283045989816821, -0.6711713557877063], [-2.0169579203546326, 3.524455258871704, -0.6214111693177178], [0.9809271994195118, 2.1566775839163737, -1.7135621949231843], [-0.9145558447621661, 2.4557522019201006, -2.72059036682588], [-2.0119541879597143, 0.6042709716797757, -1.0485616399676005], [-0.4099470658637756, -0.11587296658231545, -2.3181096973944277], [1.5770424436171668, 0.0, 0.0], [2.2927181468939155, 1.3915527243580559, 0.0], [3.782335574419717, 1.3186147352454574, -0.4807937764724458], [4.419592497958774, 0.34925073248439187, 0.1497631852713548], [3.8473329977166397, 1.1104432852891908, -1.7795740852228743], [4.36801041105158, 2.4684340853352014, -0.20560554420871952], [2.2938967322202157, 1.8704189044736101, 1.2405689893126295], [1.6494649440008826, 2.2352148943143337, -0.8029305726284923], [1.99742245733348, -0.6906780683055226, 1.0535722235493055], [1.9277183224308931, -0.6529932317206313, -1.110224125209528], [-0.3501493572534772, 0.7132786644586384, 1.0706086973199285], [-0.42668432219275876, -1.2491488329668505, 0.15362238828850008], [1.511120578020705, 3.997671524812827, 1.331413611461606]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0292', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
