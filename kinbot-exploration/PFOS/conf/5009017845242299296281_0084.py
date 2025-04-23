import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0084'
logfile = 'conf/5009017845242299296281_0084.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863835, 0.6217394783082171, -1.2501828803164987], [-0.3466020415139052, 2.1278181305643273, -1.566386336981117], [-0.9873929842445267, 2.7656987490870764, -2.8486010727213653], [-2.462172495902625, 2.278219176447381, -3.065336434982646], [-3.096370318030975, 3.1109911304421662, -3.872049480696152], [-2.4735048074690122, 1.056833510396821, -3.5916335697209933], [-3.4143213947282436, 2.1972362419820666, -1.4502458586044864], [-4.8025540278360275, 2.1453535354178994, -1.7544651034661376], [-2.7471251593519987, 1.2585061532204413, -0.6262362860864451], [-3.069060773907278, 3.6358930191985728, -0.9102417744640505], [-0.9969507321948254, 4.088639050760534, -2.7117465508075447], [-0.2788854647228145, 2.431430603605887, -3.919447328346679], [-0.7784534903451429, 2.81577717598924, -0.5057053459652474], [0.9731582247379417, 2.2459721722410175, -1.6474522505842997], [-2.0119541879597196, 0.6042709716797822, -1.0485616399675892], [-0.4099470658637787, -0.11587296658230001, -2.318109697394428], [1.5770424436171657, 0.0, 0.0], [2.292718146893916, 1.3915527243580592, 0.0], [2.3410798567223234, 2.059852692894932, 1.4165023767064735], [1.1453916410070368, 2.033172555891593, 1.975951122830792], [3.2022179922846457, 1.4422259961598307, 2.1984030628363227], [2.7215555591494964, 3.315205977928798, 1.2745358845394033], [1.62923233919397, 2.212255867310577, -0.8090479336198855], [3.5455868300943787, 1.2600392214310816, -0.42809146886198596], [1.997422457333488, -0.6906780683055278, 1.0535722235492957], [1.9277183224308945, -0.652993231720621, -1.1102241252095315], [-0.3501493572534754, 0.7132786644586327, 1.0706086973199365], [-0.4266843221927563, -1.2491488329668516, 0.15362238828850022], [-2.2421060714615866, 3.9518184565965107, -1.3004973635385393]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0084', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
