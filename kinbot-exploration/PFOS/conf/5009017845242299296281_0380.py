import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0380'
logfile = 'conf/5009017845242299296281_0380.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586385, 0.6217394783082139, -1.250182880316502], [-0.3976197158559567, -0.07566485901595706, -2.6339101198206682], [-1.1233882121466345, 0.4940238789133676, -3.9028837901505384], [-1.223597108641377, 2.058587594546231, -3.861850684121395], [-1.448535500070833, 2.526141827921839, -5.07729204510541], [-2.2078139831012673, 2.4384512294089107, -3.0519613068047615], [0.3571901289937559, 2.842070701354042, -3.2223749130466377], [1.4594646359569514, 2.172557421464403, -3.821686090283016], [0.16526526067763223, 4.2446139743930065, -3.259711005811983], [0.24551112677953896, 2.385402646423242, -1.7192033600031542], [-0.4343595853052519, 0.14440355649341488, -4.985471769970212], [-2.353319567400158, 0.0018245094686667124, -3.9766295939055754], [0.9186702629217403, 0.05926285465343825, -2.8183083047232307], [-0.6863393184328408, -1.3666682794365146, -2.523228685784666], [-0.2560445759534321, 1.8735740976390431, -1.3916701657561947], [-2.007602477187447, 0.6427130616946731, -1.054309216628063], [1.5770424436171637, 0.0, 0.0], [2.292718146893911, 1.391552724358058, 0.0], [3.7823355744197094, 1.3186147352454698, -0.4807937764724434], [4.419592497958774, 0.34925073248440064, 0.14976318527135074], [3.8473329977166326, 1.110443285289207, -1.7795740852228783], [4.368010411051572, 2.4684340853352165, -0.20560554420871746], [2.2938967322202135, 1.8704189044736113, 1.240568989312634], [1.6494649440008757, 2.23521489431434, -0.8029305726284867], [1.9974224573334844, -0.6906780683055239, 1.0535722235492981], [1.927718322430896, -0.652993231720622, -1.1102241252095335], [-0.35014935725347646, 0.7132786644586325, 1.0706086973199342], [-0.4266843221927524, -1.249148832966852, 0.15362238828850253], [0.5916409515083425, 3.0724942942125506, -1.1326357179810407]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0380', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
