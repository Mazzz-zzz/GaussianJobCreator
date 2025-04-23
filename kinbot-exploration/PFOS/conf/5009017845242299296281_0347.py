import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0347'
logfile = 'conf/5009017845242299296281_0347.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863883, -1.3935598728845964, 0.08664925740764325], [-0.3976197158559593, -2.2432006455416404, 1.3824827499919066], [-1.1233882121466405, -3.6270084497455493, 1.5236046658601465], [-2.585410647330543, -3.463818501424388, 2.06715282268516], [-3.1611759047412464, -2.4181551748563592, 1.5000685971009284], [-3.299493051811737, -4.556084848749345, 1.8103638301937364], [-2.6135892170537485, -3.204520565863992, 3.925548334931535], [-1.577339605376771, -2.2908436562316123, 4.262909054212745], [-3.9699523510745336, -3.0381260462913007, 4.296962772148894], [-2.173316102013871, -4.651783584328496, 4.363876943816577], [-1.1808982350936326, -4.206246265904724, 0.3277084738360304], [-0.452629532866642, -4.4030364263425215, 2.3652062195649637], [0.9186702629217364, -2.470358014913697, 1.3578310147309407], [-0.6863393184328458, -1.5018460017288942, 2.445183791430715], [-0.25604457595343455, -2.1420087660532796, -0.9267276815498324], [-2.0076024771874503, -1.234415095891304, -0.029451230457638908], [1.5770424436171646, 0.0, 0.0], [2.2927181468939155, 1.3915527243580539, 0.0], [2.3410798567223208, 2.0598526928949257, 1.416502376706474], [1.145391641007035, 2.033172555891588, 1.9759511228307989], [3.202217992284648, 1.442225996159824, 2.1984030628363307], [2.7215555591494907, 3.315205977928801, 1.2745358845394037], [1.629232339193968, 2.2122558673105734, -0.8090479336198848], [3.545586830094372, 1.2600392214310834, -0.42809146886198485], [1.9974224573334838, -0.6906780683055347, 1.0535722235492913], [1.9277183224308945, -0.6529932317206288, -1.11022412520954], [-0.35014935725347524, 0.5705349971623159, -1.1530217920585786], [-0.4266843221927545, 0.757615307331302, 1.0049834283127266], [-2.6270516051456454, -4.908720267119547, 5.178781859085893]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0347', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
