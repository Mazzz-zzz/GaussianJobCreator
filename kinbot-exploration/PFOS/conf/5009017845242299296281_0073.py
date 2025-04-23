import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0073'
logfile = 'conf/5009017845242299296281_0073.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863808, 0.7718203945763855, 1.1635336229088482], [-2.2709622836291916, 0.7431123812655692, 1.1797556627388976], [-2.997023964301899, 1.4656877041954879, 2.368157397369679], [-4.475954278931121, 1.8399334351693202, 2.0044399021614083], [-4.4929097750280675, 2.9333326776743074, 1.2623532827438186], [-5.055401694144313, 0.8443639468270255, 1.3396225381902014], [-5.510548804511903, 2.1626487464021174, 3.5364752564039885], [-4.736207008264464, 2.941165239636476, 4.440210759832049], [-6.805968913875699, 2.5202637789265823, 3.089983374712476], [-5.5835413086893, 0.6851580088022946, 4.077063680669295], [-3.0150768858013444, 0.6515380066656309, 3.4197399735582446], [-2.351095681931018, 2.582573574372418, 2.677628572960432], [-2.6427825570536134, 1.3431159399560084, 0.04548369333676461], [-2.6677893778920354, -0.5232706377056658, 1.1465664052870874], [-0.371045161828279, 2.061737287215915, 1.0602591291106096], [-0.24552532002049005, 0.3056731502912632, 2.324090564665853], [1.5770424436171666, 0.0, 0.0], [2.292718146893919, 1.3915527243580545, 0.0], [2.341079856722328, 2.0598526928949243, 1.4165023767064746], [1.145391641007037, 2.0331725558915905, 1.9759511228307929], [3.2022179922846448, 1.442225996159829, 2.198403062836326], [2.7215555591495044, 3.3152059779287892, 1.2745358845394064], [1.6292323391939778, 2.212255867310578, -0.8090479336198891], [3.545586830094382, 1.2600392214310685, -0.4280914688619804], [1.9974224573334807, -0.6906780683055298, 1.0535722235492957], [1.9277183224308958, -0.6529932317206282, -1.1102241252095344], [-0.35014935725347973, -1.2838136616209417, 0.08241309473864603], [-0.4266843221927545, 0.491533525635547, -1.1586058166012259], [-6.445445029090473, 0.5182280458925361, 4.483534775505774]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0073', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
