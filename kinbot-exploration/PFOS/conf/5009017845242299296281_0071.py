import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0071'
logfile = 'conf/5009017845242299296281_0071.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863808, 0.7718203945763844, 1.1635336229088489], [-2.2709622836291916, 0.7431123812655673, 1.1797556627388992], [-2.9970239643019005, 1.465687704195485, 2.3681573973696795], [-4.4759542789311215, 1.8399334351693162, 2.0044399021614097], [-4.492909775028067, 2.933332677674306, 1.2623532827438213], [-5.055401694144313, 0.8443639468270194, 1.3396225381902023], [-5.510548804511903, 2.1626487464021142, 3.5364752564039907], [-6.681358628388656, 2.8670747476411513, 3.142242618356278], [-5.522192778809229, 0.9608503418605333, 4.285409021400533], [-4.537069407306329, 3.1665721904449393, 4.261056660751672], [-3.0150768858013444, 0.6515380066656279, 3.419739973558245], [-2.3510956819310196, 2.582573574372415, 2.6776285729604337], [-2.642782557053614, 1.3431159399560073, 0.04548369333676517], [-2.6677893778920336, -0.5232706377056688, 1.1465664052870883], [-0.3710451618282805, 2.0617372872159136, 1.0602591291106112], [-0.24552532002049005, 0.30567315029126113, 2.3240905646658536], [1.5770424436171655, 0.0, 0.0], [2.2927181468939173, 1.3915527243580577, 0.0], [2.3410798567223257, 2.0598526928949266, 1.4165023767064755], [1.1453916410070342, 2.0331725558915923, 1.9759511228307935], [3.2022179922846417, 1.4422259961598307, 2.1984030628363276], [2.7215555591495013, 3.315205977928791, 1.2745358845394077], [1.6292323391939734, 2.212255867310578, -0.8090479336198854], [3.54558683009438, 1.260039221431073, -0.4280914688619806], [1.9974224573334818, -0.6906780683055279, 1.053572223549295], [1.927718322430895, -0.6529932317206278, -1.1102241252095355], [-0.3501493572534783, -1.2838136616209417, 0.08241309473864603], [-0.4266843221927556, 0.49153352563554703, -1.1586058166012252], [-4.789453904550381, 4.079783838639699, 4.065329969040315]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0071', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
