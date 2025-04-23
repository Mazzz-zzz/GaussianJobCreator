import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0147'
logfile = 'conf/5009017845242299296281_0147.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863832, 0.6217394783082129, -1.2501828803165003], [-0.39761971585595696, -0.07566485901595933, -2.633910119820666], [-0.7491833517666592, -1.5988198371506719, -2.768004365822093], [-0.9314389667275281, -2.0235549081443134, -4.266659819649045], [-2.1269824923960106, -1.6577038200326435, -4.694831050159316], [0.003036253824343063, -1.4576751685399543, -5.025234946470934], [-0.779157606893525, -3.880764665899366, -4.488573393504993], [-1.3384194975913308, -4.225327135879142, -5.749895655936187], [0.5273541439920885, -4.238432416879575, -4.075708286775726], [-1.7766523796147589, -4.327900975178208, -3.354723346484968], [0.24079307197239405, -2.3162672121361205, -2.244356108808551], [-1.879565640212353, -1.85805637597216, -2.123440869805499], [-1.1305489784219345, 0.5952337900285993, -3.5269594944504687], [0.8878652548597775, 0.09037708117406115, -2.9207161247297275], [-0.2560445759534293, 1.8735740976390407, -1.3916701657561954], [-2.007602477187446, 0.6427130616946756, -1.0543092166280608], [1.5770424436171655, 0.0, 0.0], [2.2927181468939137, 1.3915527243580543, 0.0], [1.6005215470082481, 2.440721904563897, -0.9357086002340294], [1.3760692761371613, 1.914535025410574, -2.1257143081021437], [0.46117091336275284, 2.8519467789832476, -0.4188289776134487], [2.4048649076934905, 3.478080011182911, -1.0689303403306858], [3.540059223330463, 1.229517412484612, -0.431521055692746], [2.308846803952247, 1.896094738758376, 1.2310220414904711], [1.9974224573334856, -0.6906780683055256, 1.0535722235492981], [1.9277183224308945, -0.6529932317206273, -1.1102241252095304], [-0.35014935725347357, 0.7132786644586394, 1.0706086973199322], [-0.42668432219275637, -1.2491488329668508, 0.1536223882885035], [-2.6547110819989004, -4.493397190492406, -3.725649625656655]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0147', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
