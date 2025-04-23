import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0085'
logfile = 'conf/5009017845242299296281_0085.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863857, 0.6217394783082171, -1.2501828803165], [-0.3976197158559617, -0.07566485901595055, -2.633910119820666], [-1.1233882121466408, 0.4940238789133761, -3.9028837901505353], [-1.2235971086413793, 2.058587594546239, -3.8618506841213898], [-1.4485355000708344, 2.526141827921848, -5.0772920451054055], [-2.2078139831012695, 2.4384512294089196, -3.051961306804754], [0.3571901289937549, 2.842070701354048, -3.2223749130466324], [0.38715418254115325, 2.6928733082614187, -1.8084183008412356], [1.4000823940693825, 2.4276938864492807, -4.086101924528155], [0.022212612175698726, 4.342403123643472, -3.56491734204275], [-0.43435958530525853, 0.14440355649342373, -4.985471769970211], [-2.3533195674001637, 0.0018245094686790727, -3.9766295939055722], [0.9186702629217364, 0.05926285465344263, -2.818308304723231], [-0.6863393184328467, -1.3666682794365086, -2.5232286857846686], [-0.2560445759534322, 1.8735740976390476, -1.3916701657561903], [-2.0076024771874486, 0.6427130616946773, -1.0543092166280574], [1.5770424436171615, 0.0, 0.0], [2.292718146893912, 1.3915527243580563, 0.0], [1.6005215470082503, 2.440721904563902, -0.9357086002340247], [1.3760692761371567, 1.9145350254105815, -2.1257143081021423], [0.4611709133627504, 2.8519467789832578, -0.41882897761344184], [2.404864907693496, 3.4780800111829153, -1.0689303403306827], [3.540059223330462, 1.2295174124846104, -0.43152105569275023], [2.308846803952252, 1.8960947387583753, 1.2310220414904696], [1.9974224573334858, -0.6906780683055291, 1.0535722235492928], [1.9277183224308923, -0.6529932317206217, -1.1102241252095382], [-0.35014935725347507, 0.7132786644586334, 1.0706086973199356], [-0.4266843221927539, -1.249148832966852, 0.15362238828850028], [-0.6210906780684436, 4.386890818936618, -4.286127511427906]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0085', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
