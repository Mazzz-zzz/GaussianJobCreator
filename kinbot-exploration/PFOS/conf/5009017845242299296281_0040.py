import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0040'
logfile = 'conf/5009017845242299296281_0040.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863801, 0.7718203945763868, 1.163533622908848], [-2.2709622836291903, 0.7431123812655713, 1.1797556627388988], [-2.997023964301899, 1.4656877041954923, 2.368157397369677], [-4.475954278931118, 1.8399334351693255, 2.0044399021614088], [-4.492909775028066, 2.933332677674312, 1.2623532827438178], [-5.0554016941443125, 0.8443639468270289, 1.3396225381902016], [-5.510548804511902, 2.1626487464021205, 3.5364752564039863], [-5.849703873439191, 0.9091650695481511, 4.116288604725304], [-4.883527770796265, 3.218954161488872, 4.241046235784067], [-6.787288074297166, 2.7535783494098443, 2.82846277693323], [-3.015076885801341, 0.6515380066656374, 3.4197399735582454], [-2.351095681931016, 2.582573574372421, 2.6776285729604306], [-2.642782557053612, 1.3431159399560106, 0.04548369333676349], [-2.6677893778920327, -0.5232706377056646, 1.1465664052870899], [-0.371045161828278, 2.0617372872159163, 1.060259129110609], [-0.24552532002049057, 0.3056731502912632, 2.324090564665853], [1.5770424436171668, 0.0, 0.0], [2.2927181468939204, 1.3915527243580548, 0.0], [2.3410798567223323, 2.0598526928949257, 1.4165023767064733], [1.145391641007042, 2.0331725558915945, 1.975951122830791], [3.2022179922846465, 1.4422259961598276, 2.1984030628363254], [2.7215555591495146, 3.31520597792879, 1.274535884539404], [1.6292323391939811, 2.212255867310578, -0.8090479336198891], [3.5455868300943854, 1.2600392214310667, -0.4280914688619814], [1.9974224573334813, -0.6906780683055296, 1.0535722235492957], [1.927718322430895, -0.6529932317206288, -1.1102241252095344], [-0.35014935725348006, -1.2838136616209415, 0.08241309473864715], [-0.42668432219275354, 0.4915335256355479, -1.158605816601226], [-7.448689843050857, 2.0602682562226584, 2.694929254287201]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0040', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
