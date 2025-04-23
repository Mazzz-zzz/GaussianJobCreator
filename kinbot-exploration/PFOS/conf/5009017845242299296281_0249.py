import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0249'
logfile = 'conf/5009017845242299296281_0249.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863834, -1.3935598728845984, 0.08664925740765084], [-0.3466020415139029, -2.4204394252486643, -1.059551387211261], [1.1624535490467058, -2.7999021924931093, -1.2603960465750452], [1.3229555267971613, -4.17169637078313, -2.0033664407024308], [2.5379668391891173, -4.265983898549034, -2.5145995733771147], [1.1168882370515762, -5.1808849811579964, -1.1619541390693033], [0.08826989713528051, -4.35440177328728, -3.404734865820267], [0.5381926157038321, -5.3970590001345595, -4.260811067118386], [-1.2025433687641451, -4.305859813711145, -2.8243887945950163], [0.34402262312628656, -2.972098158293413, -4.1149249030159885], [1.750051457246916, -1.8550381832545337, -1.989060753302593], [1.7633282149005842, -2.8945149144395588, -0.08123199086868525], [-1.0117971004258974, -3.5349698356523596, -0.7430628651027303], [-0.821528296712381, -1.9466956411138778, -2.2051128650057987], [-2.0119541879597174, -1.2102165034856969, 0.0009668077395966845], [-0.4099470658637775, -1.9496054034114845, 1.2594037813693517], [1.577042443617167, 0.0, 0.0], [2.2927181468939164, 1.391552724358056, 0.0], [3.782335574419717, 1.3186147352454594, -0.48079377647244387], [4.41959249795878, 0.3492507324843892, 0.1497631852713513], [3.84733299771664, 1.1104432852891992, -1.779574085222876], [4.368010411051582, 2.4684340853352, -0.20560554420871474], [2.29389673222022, 1.8704189044736053, 1.2405689893126328], [1.6494649440008833, 2.235214894314338, -0.80293057262848], [1.9974224573334833, -0.6906780683055272, 1.0535722235492984], [1.927718322430895, -0.6529932317206247, -1.1102241252095315], [-0.35014935725347296, 0.5705349971623099, -1.1530217920585821], [-0.4266843221927547, 0.7576153073313033, 1.0049834283127284], [0.7028940834630952, -2.3318004857815438, -3.484692761477961]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0249', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
