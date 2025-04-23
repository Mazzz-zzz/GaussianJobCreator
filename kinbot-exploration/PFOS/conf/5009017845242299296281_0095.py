import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0095'
logfile = 'conf/5009017845242299296281_0095.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863819, 0.7718203945763858, 1.1635336229088487], [-2.270962283629193, 0.743112381265566, 1.1797556627388996], [-3.020318448930584, 1.3845311280592587, -0.040337232539326216], [-4.480058160568045, 0.8292457901518767, -0.18317744667967537], [-4.45376171327607, -0.3598415970703149, -0.7592691281730983], [-5.056700168542038, 0.7301262577332812, 1.0112091546709316], [-5.555811090753356, 1.9566623859997165, -1.2288116005401828], [-5.935868144708009, 3.0723663321302475, -0.4330777580161434], [-4.931749541331282, 2.0615750382471933, -2.495657009254022], [-6.7987732354635675, 1.0024547339658836, -1.3871640347698775], [-3.086116569243985, 2.7007512384555907, 0.1391809583425367], [-2.3639288743924594, 1.1178871860085962, -1.1621587584935829], [-2.5961406643712652, -0.5518249635012422, 1.2269126406770388], [-2.6897770261666514, 1.332584330848514, 2.293046394029769], [-0.37104516182828207, 2.0617372872159105, 1.0602591291106105], [-0.24552532002049077, 0.30567315029126346, 2.324090564665855], [1.5770424436171673, 0.0, 0.0], [2.2927181468939173, 1.3915527243580528, 0.0], [3.7823355744197165, 1.3186147352454567, -0.4807937764724417], [4.419592497958776, 0.34925073248439076, 0.14976318527134935], [3.8473329977166415, 1.1104432852891943, -1.7795740852228759], [4.368010411051579, 2.4684340853352067, -0.20560554420871519], [2.2938967322202157, 1.8704189044736046, 1.2405689893126346], [1.6494649440008815, 2.235214894314333, -0.8029305726284865], [1.9974224573334847, -0.6906780683055256, 1.053572223549306], [1.9277183224308967, -0.6529932317206291, -1.1102241252095264], [-0.3501493572534747, -1.283813661620943, 0.08241309473865067], [-0.42668432219275393, 0.4915335256355453, -1.1586058166012254], [-7.468324495646526, 1.2092449236636011, -0.7201548038530967]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0095', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
