import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0220'
logfile = 'conf/5009017845242299296281_0220.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863829, 0.7718203945763793, 1.163533622908848], [-0.3466020415139028, 0.29262129468434017, 2.6259377241923785], [1.1624535490467065, 0.30841610108309514, 3.054984450098302], [1.3229555267971613, 0.3508819546540548, 4.6144782543247524], [2.5379668391891173, -0.044715161615564956, 4.951750214967398], [1.1168882370515751, 1.5841606881125097, 5.0677550773027376], [0.08826989713528173, -0.7713860003073056, 5.4733899868609255], [0.026477469084687354, -1.9970674818391692, 4.754837546550922], [0.345434384503417, -0.6784482742478193, 6.8628451956242715], [-1.217944947265287, 0.05603896078061508, 5.173868128150664], [1.7500514572469168, -0.7950580504033847, 2.601040568339856], [1.7633282149005822, 1.3769084895275159, 2.547339442971939], [-1.0117971004258977, 1.1239736000383669, 3.432905111838008], [-0.8215282967123797, -0.9363359387499672, 2.788444311143953], [-2.0119541879597165, 0.6059455318059099, 1.0475948322279918], [-0.40994706586377916, 2.065478369993783, 1.0587059160250758], [1.5770424436171682, 0.0, 0.0], [2.2927181468939137, 1.3915527243580603, 0.0], [2.34107985672232, 2.0598526928949337, 1.4165023767064715], [1.1453916410070288, 2.0331725558915905, 1.9759511228307938], [3.2022179922846385, 1.4422259961598354, 2.1984030628363307], [2.7215555591494875, 3.315205977928801, 1.2745358845394068], [1.62923233919397, 2.2122558673105783, -0.8090479336198886], [3.5455868300943765, 1.2600392214310878, -0.42809146886197835], [1.9974224573334873, -0.6906780683055237, 1.0535722235492986], [1.927718322430897, -0.6529932317206256, -1.1102241252095335], [-0.350149357253473, -1.2838136616209446, 0.08241309473864965], [-0.42668432219275404, 0.4915335256355439, -1.1586058166012294], [-1.8300838645810733, -0.00021205134213187835, 5.920918431837521]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0220', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
