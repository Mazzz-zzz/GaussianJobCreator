import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0331'
logfile = 'conf/5009017845242299296281_0331.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863824, 0.6217394783082101, -1.250182880316507], [-0.34660204151389995, 2.1278181305643167, -1.5663863369811326], [-0.7363023803695176, 3.202618946498, -0.4919021834050533], [-0.5090926255515525, 2.674725836015372, 0.9673064206754405], [0.6213430942919618, 1.9927921769719528, 1.0271688549938693], [-0.45968485817329485, 3.6902998099574162, 1.8246677310500485], [-1.8990683064763108, 1.5434550264286835, 1.523954687473987], [-2.220326922430838, 0.6708970121473475, 0.44791267100883553], [-1.574188825241977, 1.105243334171232, 2.830775488919148], [-3.0356793420946144, 2.6260608394842855, 1.652883493634265], [0.018370741190135768, 4.283045989816821, -0.6711713557877063], [-2.0169579203546317, 3.5244552588717033, -0.6214111693177192], [0.9809271994195126, 2.1566775839163705, -1.713562194923185], [-0.9145558447621657, 2.455752201920098, -2.72059036682588], [-2.0119541879597134, 0.604270971679775, -1.0485616399676003], [-0.4099470658637681, -0.11587296658231455, -2.3181096973944304], [1.5770424436171675, 0.0, 0.0], [2.2927181468939155, 1.3915527243580548, 0.0], [2.341079856722321, 2.0598526928949266, 1.4165023767064786], [1.1453916410070335, 2.0331725558915945, 1.9759511228307913], [3.2022179922846385, 1.4422259961598245, 2.19840306283633], [2.7215555591495004, 3.3152059779287937, 1.274535884539408], [1.629232339193981, 2.2122558673105797, -0.8090479336198884], [3.545586830094382, 1.2600392214310743, -0.42809146886197696], [1.9974224573334822, -0.690678068305522, 1.0535722235493061], [1.927718322430896, -0.65299323172063, -1.1102241252095266], [-0.3501493572534767, 0.7132786644586381, 1.070608697319928], [-0.42668432219275493, -1.2491488329668512, 0.15362238828850242], [-3.617346185848723, 2.4156852018916477, 2.3967706898412073]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0331', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
