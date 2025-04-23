import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0353'
logfile = 'conf/5009017845242299296281_0353.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863824, 0.6217394783082101, -1.250182880316507], [-0.34660204151389995, 2.1278181305643167, -1.5663863369811326], [-0.7363023803695176, 3.202618946498, -0.4919021834050533], [-0.5090926255515525, 2.674725836015372, 0.9673064206754405], [0.6213430942919618, 1.9927921769719528, 1.0271688549938693], [-0.45968485817329485, 3.6902998099574162, 1.8246677310500485], [-1.8990683064763108, 1.5434550264286835, 1.523954687473987], [-2.220326922430838, 0.6708970121473475, 0.44791267100883553], [-1.574188825241977, 1.105243334171232, 2.830775488919148], [-3.0356793420946144, 2.6260608394842855, 1.652883493634265], [0.018370741190135768, 4.283045989816821, -0.6711713557877063], [-2.0169579203546317, 3.5244552588717033, -0.6214111693177192], [0.9809271994195126, 2.1566775839163705, -1.713562194923185], [-0.9145558447621657, 2.455752201920098, -2.72059036682588], [-2.0119541879597134, 0.604270971679775, -1.0485616399676003], [-0.4099470658637681, -0.11587296658231455, -2.3181096973944304], [1.5770424436171675, 0.0, 0.0], [2.2927181468939155, 1.3915527243580548, 0.0], [1.600521547008245, 2.4407219045638975, -0.9357086002340314], [1.3760692761371671, 1.9145350254105737, -2.125714308102146], [0.4611709133627493, 2.851946778983244, -0.41882897761345483], [2.404864907693491, 3.478080011182911, -1.068930340330686], [3.5400592233304677, 1.2295174124846144, -0.43152105569274196], [2.3088468039522447, 1.896094738758375, 1.231022041490469], [1.9974224573334824, -0.690678068305522, 1.053572223549306], [1.9277183224308962, -0.6529932317206301, -1.1102241252095266], [-0.3501493572534767, 0.7132786644586381, 1.070608697319928], [-0.42668432219275493, -1.2491488329668512, 0.15362238828850242], [-3.5663728929085403, 2.654067331790809, 0.8444638141303626]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0353', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
