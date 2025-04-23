import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0235'
logfile = 'conf/5009017845242299296281_0235.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863807, 0.7718203945763862, 1.1635336229088493], [-2.270962283629192, 0.7431123812655699, 1.1797556627389014], [-2.9699917885127123, -0.6604598778102058, 1.233811274632789], [-4.434530032795169, -0.559565237651802, 1.7856565080301303], [-5.030583914311105, 0.5082304477086564, 1.2846234874293212], [-5.127712885957739, -1.6471563109736436, 1.4607823325943794], [-4.4665010348623895, -0.41736338044220905, 3.656597399118879], [-4.226968110733983, -1.7075006466777913, 4.2048681537088655], [-3.7215317568045045, 0.7387617773290985, 3.993869049226321], [-5.992464116672816, -0.06933155880081121, 3.832248057772909], [-3.0171454311974584, -1.164914771591805, 0.0040579064213605645], [-2.2839557356380094, -1.4745107951449294, 2.0257111670633865], [-2.6212054717929205, 1.418984801741657, 2.2776416618875865], [-2.713377980741126, 1.397124924610703, 0.11259346120335914], [-0.37104516182827757, 2.061737287215914, 1.0602591291106138], [-0.24552532002048638, 0.3056731502912653, 2.324090564665853], [1.577042443617164, 0.0, 0.0], [2.292718146893918, 1.3915527243580552, 0.0], [1.6005215470082532, 2.4407219045638975, -0.935708600234027], [1.3760692761371622, 1.9145350254105797, -2.1257143081021415], [0.4611709133627553, 2.851946778983256, -0.41882897761344673], [2.4048649076935016, 3.47808001118291, -1.0689303403306831], [3.540059223330462, 1.229517412484601, -0.43152105569274757], [2.308846803952253, 1.8960947387583698, 1.2310220414904673], [1.997422457333483, -0.6906780683055318, 1.053572223549296], [1.9277183224308907, -0.6529932317206231, -1.11022412520954], [-0.35014935725347984, -1.2838136616209443, 0.08241309473864619], [-0.42668432219275787, 0.4915335256355502, -1.1586058166012263], [-6.34051051828528, 0.33320452941809425, 3.024293364301768]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0235', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
