import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0120'
logfile = 'conf/5009017845242299296281_0120.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863837, 0.621739478308217, -1.250182880316503], [-0.3976197158559589, -0.07566485901595643, -2.633910119820667], [-0.749183351766663, -1.5988198371506686, -2.7680043658220956], [-0.9314389667275326, -2.0235549081443094, -4.266659819649047], [-0.8151262318136152, -3.3351810806509357, -4.378954407094112], [-2.1256098696188777, -1.6449135451208223, -4.713299062717092], [0.35577392705740013, -1.2298774894092341, -5.377885472781138], [0.3870193333205283, -1.9474277222107508, -6.60531221022983], [0.17369960609959667, 0.17122805603693517, -5.2824445810431175], [1.6380449215560549, -1.6020433096845461, -4.542531458787524], [0.24079307197238695, -2.31626721213612, -2.2443561088085526], [-1.8795656402123568, -1.8580563759721536, -2.1234408698055014], [-1.130548978421933, 0.5952337900286033, -3.5269594944504687], [0.8878652548597777, 0.09037708117406118, -2.9207161247297284], [-0.25604457595342817, 1.873574097639044, -1.3916701657561978], [-2.0076024771874463, 0.642713061694682, -1.054309216628064], [1.5770424436171622, 0.0, 0.0], [2.2927181468939164, 1.3915527243580512, 0.0], [1.6005215470082588, 2.4407219045638975, -0.935708600234027], [1.3760692761371698, 1.9145350254105806, -2.125714308102142], [0.46117091336276284, 2.8519467789832555, -0.4188289776134475], [2.4048649076935043, 3.478080011182908, -1.0689303403306845], [3.5400592233304664, 1.2295174124846024, -0.431521055692743], [2.3088468039522523, 1.8960947387583693, 1.2310220414904731], [1.9974224573334818, -0.6906780683055276, 1.0535722235492986], [1.927718322430891, -0.6529932317206288, -1.1102241252095308], [-0.35014935725347246, 0.7132786644586406, 1.0706086973199318], [-0.4266843221927601, -1.2491488329668479, 0.153622388288502], [2.290171943046478, -0.8891534426249367, -4.59233880745981]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0120', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
