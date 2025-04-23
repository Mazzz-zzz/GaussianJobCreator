import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0009'
logfile = 'conf/5009017845242299296281_0009.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863835, -1.3935598728845968, 0.08664925740765074], [-0.3466020415139052, -2.420439425248664, -1.0595513872112619], [-0.9873929842445267, -3.849810268767844, -0.9708648395635668], [-2.462172495902625, -3.793768812064689, -0.44032746470096895], [-3.0963703180309756, -4.908788780214294, -0.758172609562905], [-2.4735048074690122, -3.638862667661778, 0.8805721172861672], [-3.4143213947282436, -2.354567876275694, -1.1777394743700782], [-4.8025540278360275, -2.592088117363921, -0.980698110037591], [-2.7471251593519987, -1.1715896091326994, -0.7767801564647083], [-3.0690607739072786, -2.6062390098709813, -2.6936548328364385], [-0.9969507321948254, -4.392760927004431, -2.1849920094599384], [-0.2788854647228145, -4.610056256946216, -0.145957006088289], [-0.7784534903451429, -1.8458422644301251, -2.185681892820462], [0.9731582247379417, -2.54972158664836, -1.1213428320614889], [-2.0119541879597196, -1.2102165034856958, 0.0009668077395966837], [-0.4099470658637787, -1.9496054034114814, 1.2594037813693524], [1.5770424436171657, 0.0, 0.0], [2.2927181468939186, 1.391552724358055, 0.0], [3.7823355744197125, 1.3186147352454622, -0.4807937764724419], [4.4195924979587735, 0.34925073248439276, 0.14976318527134758], [3.8473329977166366, 1.110443285289201, -1.77957408522288], [4.36801041105158, 2.468434085335202, -0.2056055442087168], [2.2938967322202175, 1.870418904473606, 1.2405689893126326], [1.6494649440008822, 2.235214894314338, -0.8029305726284864], [1.9974224573334842, -0.6906780683055272, 1.053572223549297], [1.9277183224308936, -0.6529932317206243, -1.110224125209535], [-0.35014935725347596, 0.5705349971623112, -1.1530217920585806], [-0.4266843221927563, 0.7576153073313034, 1.0049834283127286], [-3.7745605871978203, -3.1132335616509716, -3.1193177565312493]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0009', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
