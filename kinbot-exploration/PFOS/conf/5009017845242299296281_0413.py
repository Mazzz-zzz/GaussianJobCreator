import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0413'
logfile = 'conf/5009017845242299296281_0413.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863852, 0.6217394783082121, -1.2501828803165], [-2.2709622836291934, 0.6501421835576475, -1.233432031412176], [-2.9970239643019028, 1.3180406141844365, -2.45340148453262], [-2.2704853996782592, 1.0013618477382582, -3.8067090741066765], [-1.213602341516464, 1.7830531864924313, -3.942572223504435], [-1.8869755643483324, -0.27172278628530283, -3.8395078754768472], [-3.3935951456019624, 1.282594210680647, -5.283596404866892], [-4.105849023147379, 2.495878636846642, -5.076072134106644], [-2.628917663197577, 0.9996625578462482, -6.441428344042906], [-4.377248131761582, 0.07403875019459004, -5.054766486461104], [-4.239914773111563, 0.8496739297323483, -2.5230192895671477], [-3.023054582933855, 2.635152701311813, -2.2953889586822838], [-2.6427825570536156, -0.6321679360904267, -1.1859143708980975], [-2.667789377892034, 1.254590952957258, -0.12011753733594914], [-0.371045161828282, -0.11265730320380403, -2.3156464312139], [-0.24552532002049046, 1.8598848945507154, -1.4267659957399805], [1.5770424436171642, 0.0, 0.0], [2.292718146893916, 1.391552724358059, 0.0], [2.34107985672232, 2.059852692894926, 1.4165023767064744], [1.1453916410070384, 2.033172555891597, 1.9759511228307964], [3.202217992284648, 1.442225996159831, 2.1984030628363245], [2.721555559149498, 3.3152059779287955, 1.2745358845394037], [1.6292323391939711, 2.2122558673105748, -0.809047933619884], [3.5455868300943774, 1.2600392214310723, -0.42809146886198157], [1.9974224573334807, -0.6906780683055259, 1.053572223549303], [1.927718322430895, -0.6529932317206245, -1.1102241252095284], [-0.35014935725347435, 0.7132786644586341, 1.0706086973199322], [-0.42668432219275576, -1.2491488329668528, 0.1536223882885026], [-3.9423703820558287, -0.6204849914423879, -4.540489808194015]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0413', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
