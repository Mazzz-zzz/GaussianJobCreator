import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0476'
logfile = 'conf/5009017845242299296281_0476.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586382, 0.7718203945763829, 1.163533622908851], [-2.2709622836291916, 0.7431123812655621, 1.1797556627389025], [-3.0203184489305848, 1.3845311280592574, -0.04033723253932064], [-4.4800581605680465, 0.82924579015187, -0.18317744667966848], [-5.202953426281216, 1.651766877071618, -0.9230311075200033], [-4.464154537777425, -0.3753605633517404, -0.7467115170195937], [-5.325329553432012, 0.6531877726236972, 1.4830107713752003], [-5.0309972643683745, 1.8103381033937131, 2.255540785952889], [-6.635959680959697, 0.17970425591912034, 1.2311798759833619], [-4.481042378252227, -0.5540860981190835, 2.039977194651937], [-3.086116569243987, 2.700751238455587, 0.1391809583425437], [-2.3639288743924594, 1.1178871860085953, -1.1621587584935797], [-2.5961406643712626, -0.551824963501244, 1.2269126406770399], [-2.6897770261666496, 1.3325843308485081, 2.2930463940297736], [-0.3710451618282817, 2.061737287215911, 1.060259129110611], [-0.2455253200204883, 0.3056731502912637, 2.3240905646658567], [1.5770424436171655, 0.0, 0.0], [2.2927181468939164, 1.3915527243580526, 0.0], [1.6005215470082548, 2.440721904563895, -0.935708600234029], [1.3760692761371702, 1.9145350254105822, -2.1257143081021415], [0.4611709133627595, 2.8519467789832547, -0.4188289776134507], [2.404864907693508, 3.478080011182909, -1.0689303403306805], [3.540059223330469, 1.2295174124846038, -0.43152105569274435], [2.308846803952251, 1.8960947387583702, 1.2310220414904738], [1.997422457333486, -0.6906780683055234, 1.0535722235493008], [1.9277183224308974, -0.6529932317206238, -1.110224125209531], [-0.35014935725347457, -1.283813661620944, 0.08241309473865073], [-0.4266843221927561, 0.4915335256355469, -1.158605816601225], [-3.7575475330199386, -0.23056193861479377, 2.5948192567158674]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0476', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
