import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0097'
logfile = 'conf/5009017845242299296281_0097.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586382, 0.7718203945763852, 1.1635336229088478], [-2.270962283629193, 0.7431123812655701, 1.1797556627388968], [-2.969991788512713, -0.6604598778102067, 1.2338112746327854], [-3.0567398040787164, -1.3216748454235765, -0.18564450396738305], [-4.057927071959132, -0.792630756448845, -0.866994890796395], [-1.9208515428855046, -1.1373859158669182, -0.8524994275297249], [-3.343878254277271, -3.172954158722334, -0.07619247938738649], [-3.8032649569633916, -3.6274572122401456, -1.343008109046317], [-2.232158454102362, -3.7217198011322314, 0.6081262528255562], [-4.559261845011894, -3.155833489648285, 0.9253524471498481], [-2.2659422418234927, -1.4610539209006586, 2.029018773642186], [-4.2028972120978745, -0.5310943056708308, 1.70671566954612], [-2.6212054717929245, 1.4189848017416598, 2.2776416618875777], [-2.713377980741126, 1.3971249246107025, 0.11259346120335348], [-0.3710451618282794, 2.0617372872159145, 1.0602591291106083], [-0.24552532002049035, 0.305673150291265, 2.324090564665851], [1.577042443617164, 0.0, 0.0], [2.2927181468939195, 1.3915527243580508, 0.0], [2.3410798567223274, 2.059852692894926, 1.4165023767064713], [1.1453916410070484, 2.0331725558915976, 1.975951122830795], [3.202217992284648, 1.4422259961598258, 2.198403062836327], [2.7215555591495164, 3.315205977928792, 1.2745358845393986], [1.6292323391939827, 2.212255867310574, -0.8090479336198897], [3.5455868300943862, 1.2600392214310625, -0.42809146886198696], [1.9974224573334818, -0.6906780683055308, 1.0535722235492957], [1.9277183224308942, -0.6529932317206303, -1.1102241252095313], [-0.35014935725347784, -1.2838136616209423, 0.08241309473865063], [-0.42668432219275526, 0.49153352563554686, -1.1586058166012283], [-5.395287609971865, -3.1931254597200427, 0.439932406373551]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0097', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
