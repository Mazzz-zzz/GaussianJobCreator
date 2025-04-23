import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0090'
logfile = 'conf/5009017845242299296281_0090.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863854, 0.6217394783082142, -1.2501828803165027], [-0.3976197158559566, -0.07566485901595763, -2.633910119820668], [-0.749183351766658, -1.5988198371506697, -2.768004365822095], [-2.081987459715323, -1.9555730430184126, -2.0223868238206655], [-2.9765439628333943, -1.0033833036566053, -2.2213976427734505], [-2.563171468208477, -3.1145119714526817, -2.4630936999675], [-1.8210733507287753, -2.1287237062650552, -0.17208792724638705], [-1.2114586234507911, -3.3895648474906657, 0.07505678479228862], [-1.3031811287452222, -0.8934391866980015, 0.2874951161472201], [-3.332388801547781, -2.204189361974298, 0.26468795662671124], [-0.8953137364216069, -1.8945036554923973, -4.0564951662384], [0.22590414679516438, -2.333069003291661, -2.2478963639045046], [-1.1305489784219311, 0.5952337900286009, -3.5269594944504687], [0.8878652548597795, 0.09037708117406244, -2.9207161247297275], [-0.25604457595343105, 1.8735740976390431, -1.3916701657561947], [-2.00760247718745, 0.6427130616946739, -1.0543092166280643], [1.5770424436171646, 0.0, 0.0], [2.292718146893914, 1.3915527243580577, 0.0], [2.3410798567223217, 2.0598526928949297, 1.4165023767064748], [1.145391641007032, 2.0331725558915936, 1.9759511228307902], [3.2022179922846448, 1.4422259961598287, 2.1984030628363245], [2.721555559149491, 3.3152059779287972, 1.2745358845394066], [1.6292323391939767, 2.2122558673105805, -0.8090479336198868], [3.54558683009438, 1.2600392214310787, -0.4280914688619807], [1.9974224573334838, -0.6906780683055244, 1.0535722235492981], [1.9277183224308947, -0.6529932317206257, -1.1102241252095308], [-0.35014935725347723, 0.7132786644586374, 1.0706086973199334], [-0.4266843221927588, -1.2491488329668499, 0.15362238828850225], [-3.8949003173690895, -1.776316178643826, -0.3959690801564201]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0090', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
