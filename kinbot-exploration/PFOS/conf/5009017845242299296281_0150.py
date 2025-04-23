import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0150'
logfile = 'conf/5009017845242299296281_0150.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863839, 0.7718203945763815, 1.163533622908849], [-2.270962283629191, 0.7431123812655616, 1.1797556627389043], [-2.9970239643018997, 1.4656877041954797, 2.368157397369685], [-3.0516123514917344, 0.5642090685426583, 3.650320232756651], [-1.9012352660876428, -0.06783254007289861, 3.8041315113835754], [-3.2973804088622845, 1.305110749334826, 4.727150210946642], [-4.411157127461302, -0.7245766174249222, 3.5392671423272173], [-5.653942072382576, -0.08348488425822553, 3.7979522368965464], [-4.145375387517529, -1.5044078441720634, 2.387494553596906], [-4.014070723755674, -1.5832180667783868, 4.798473193477429], [-2.3329497198238167, 2.5793667469040993, 2.6643261670914047], [-4.241488402883601, 1.7657685081001782, 2.019317070396958], [-2.6427825570536174, 1.3431159399560015, 0.04548369333676916], [-2.6677893778920323, -0.5232706377056753, 1.1465664052870943], [-0.3710451618282831, 2.0617372872159128, 1.060259129110613], [-0.2455253200204879, 0.30567315029126113, 2.3240905646658536], [1.5770424436171655, 0.0, 0.0], [2.2927181468939155, 1.3915527243580577, 0.0], [3.782335574419711, 1.3186147352454614, -0.480793776472449], [4.419592497958776, 0.3492507324843972, 0.14976318527134774], [3.8473329977166357, 1.1104432852892017, -1.7795740852228776], [4.368010411051574, 2.4684340853352102, -0.20560554420872507], [2.293896732220214, 1.8704189044736101, 1.240568989312627], [1.6494649440008777, 2.235214894314336, -0.8029305726284898], [1.9974224573334842, -0.6906780683055271, 1.0535722235492957], [1.9277183224308958, -0.6529932317206266, -1.110224125209533], [-0.35014935725347535, -1.2838136616209457, 0.08241309473864857], [-0.426684322192759, 0.4915335256355461, -1.1586058166012265], [-4.186057659751961, -2.520884017947233, 4.633653763308452]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0150', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
