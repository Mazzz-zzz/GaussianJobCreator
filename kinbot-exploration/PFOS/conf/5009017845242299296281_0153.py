import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0153'
logfile = 'conf/5009017845242299296281_0153.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863854, -1.3935598728845993, 0.08664925740764841], [-0.3976197158559566, -2.24320064554164, 1.382482749991912], [-0.749183351766658, -1.5977521800128425, 2.768620777958022], [-2.081987459715323, -0.7736518441984246, 2.704769346120314], [-2.9765439628333943, -1.4220951387203828, 1.9796542520864955], [-2.563171468208477, -0.5758457303469363, 3.9287933376525244], [-1.8210733507287753, 0.9153293364525414, 1.9295727708868982], [-1.2114586234507911, 1.7597835061018243, 2.8979208733055066], [-1.3031811287452217, 0.695697667396449, 0.6299934743433695], [-3.332388801547781, 1.3313211755016672, 1.7765400039078056], [-0.8953137364216069, -2.5657760365450497, 3.6689358763380886], [0.22590414679516438, -0.780200854570151, 3.144445207584868], [-1.1305489784219311, -3.3520534153271324, 1.2479921638695595], [0.8878652548597796, -2.5746029018458203, 1.3820892141482264], [-0.25604457595343105, -2.142008766053283, -0.9267276815498285], [-2.00760247718745, -1.2344150958913116, -0.029451230457635244], [1.5770424436171646, 0.0, 0.0], [2.292718146893912, 1.3915527243580565, 0.0], [2.341079856722316, 2.0598526928949292, 1.4165023767064748], [1.1453916410070282, 2.0331725558915927, 1.9759511228307944], [3.2022179922846403, 1.4422259961598338, 2.198403062836327], [2.7215555591494858, 3.3152059779288017, 1.274535884539406], [1.6292323391939618, 2.212255867310577, -0.8090479336198818], [3.5455868300943725, 1.2600392214310856, -0.42809146886198546], [1.9974224573334824, -0.690678068305528, 1.0535722235492966], [1.9277183224308927, -0.6529932317206235, -1.1102241252095368], [-0.3501493572534773, 0.5705349971623132, -1.1530217920585808], [-0.4266843221927588, 0.7576153073313016, 1.0049834283127277], [-3.6101933920399523, 1.8714150440715573, 2.5295908715991513]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0153', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
