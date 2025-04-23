import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0013'
logfile = 'conf/5009017845242299296281_0013.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863824, 0.6217394783082114, -1.250182880316503], [-0.34660204151390517, 2.1278181305643216, -1.5663863369811248], [-0.7363023803695233, 3.2026189464980006, -0.49190218340504316], [0.11232051632788764, 4.51250111182252, -0.6456326331137674], [0.2504583136006516, 4.810364670493017, -1.925748718653581], [-0.47946595887502164, 5.524080936392357, -0.016828158858708243], [1.8265739031510881, 4.323642832150676, 0.09418572775568784], [1.7145076756829611, 4.425991285687421, 1.5081860681126809], [2.446769046148746, 3.241378914280436, -0.5762018976717815], [2.461283180520462, 5.65436513358516, -0.4597681687873024], [-2.0209682759934693, 3.513213441320752, -0.640751600308594], [-0.5312999584009233, 2.716571178249477, 0.7255383184087576], [0.9809271994195065, 2.1566775839163776, -1.7135621949231796], [-0.9145558447621711, 2.455752201920104, -2.720590366825872], [-2.0119541879597174, 0.6042709716797755, -1.048561639967597], [-0.4099470658637744, -0.11587296658230682, -2.31810969739443], [1.577042443617164, 0.0, 0.0], [2.292718146893912, 1.3915527243580597, 0.0], [2.3410798567223177, 2.059852692894931, 1.4165023767064775], [1.1453916410070284, 2.0331725558915936, 1.9759511228307909], [3.2022179922846394, 1.442225996159833, 2.1984030628363285], [2.7215555591494898, 3.3152059779288017, 1.274535884539406], [1.629232339193969, 2.2122558673105797, -0.8090479336198877], [3.5455868300943756, 1.2600392214310825, -0.42809146886198013], [1.9974224573334827, -0.6906780683055218, 1.0535722235493008], [1.9277183224308965, -0.6529932317206227, -1.11022412520953], [-0.35014935725347873, 0.7132786644586351, 1.0706086973199318], [-0.4266843221927586, -1.249148832966852, 0.15362238828850028], [1.9952003216854641, 5.9425077211713875, -1.2570780236212342]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0013', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
