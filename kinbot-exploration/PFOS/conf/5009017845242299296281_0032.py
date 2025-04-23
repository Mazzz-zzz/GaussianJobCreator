import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0032'
logfile = 'conf/5009017845242299296281_0032.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863839, -1.3935598728845964, 0.08664925740764823], [-0.3976197158559582, -2.243200645541639, 1.3824827499919112], [-1.123388212146639, -3.627008449745548, 1.5236046658601499], [-2.5854106473305416, -3.463818501424388, 2.0671528226851645], [-2.5604417585430426, -3.2859249873038996, 3.3764398024201547], [-3.1805547950289808, -2.4237675353049792, 1.490093491344449], [-3.6360803989189696, -4.976515107083823, 1.707284258569699], [-2.8559764250854114, -6.1344368241106855, 1.9776286593932508], [-4.91376663469693, -4.732305001292866, 2.266959234588026], [-3.758030186975003, -4.802764753073875, 0.14668126822157335], [-1.1808982350936252, -4.206246265904727, 0.32770847383603807], [-0.4526295328666376, -4.403036426342518, 2.3652062195649717], [0.9186702629217389, -2.470358014913695, 1.3578310147309456], [-0.6863393184328441, -1.5018460017288948, 2.445183791430716], [-0.256044575953429, -2.142008766053283, -0.9267276815498295], [-2.0076024771874477, -1.2344150958913078, -0.02945123045764119], [1.577042443617165, 0.0, 0.0], [2.292718146893918, 1.391552724358054, 0.0], [3.782335574419717, 1.3186147352454607, -0.4807937764724353], [4.419592497958776, 0.34925073248438954, 0.14976318527135307], [3.847332997716641, 1.1104432852891988, -1.7795740852228759], [4.368010411051582, 2.468434085335201, -0.20560554420871308], [2.2938967322202157, 1.870418904473602, 1.240568989312635], [1.6494649440008855, 2.2352148943143413, -0.8029305726284788], [1.9974224573334831, -0.6906780683055311, 1.0535722235492972], [1.9277183224308965, -0.6529932317206258, -1.1102241252095342], [-0.35014935725347374, 0.5705349971623116, -1.1530217920585812], [-0.42668432219275637, 0.7576153073313003, 1.0049834283127244], [-3.0852083056366904, -5.332082781453226, -0.30397994324550187]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0032', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
