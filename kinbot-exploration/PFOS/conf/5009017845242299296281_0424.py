import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0424'
logfile = 'conf/5009017845242299296281_0424.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586381, 0.6217394783082208, -1.2501828803164983], [-0.3466020415139008, 2.127818130564333, -1.5663863369811128], [-0.9873929842445208, 2.765698749087086, -2.8486010727213604], [-0.18477695124074434, 2.3967981845477566, -4.144481096508226], [0.8897635883676692, 3.160212701690833, -4.239680081955854], [0.1812889705181365, 1.1185253607175663, -4.114242555918094], [-1.2119219398679435, 2.6451172059345125, -5.695279009542207], [-1.9175735786971912, 3.873389226204253, -5.569443336340685], [-0.382611162427967, 2.3156882622018977, -6.794859373636226], [-2.225283750583378, 1.4571193219919703, -5.489605729458032], [-2.230506898554333, 2.3112475123060654, -2.979495770641265], [-1.0036314587762258, 4.087403218990677, -2.733434203914618], [-0.7784534903451344, 2.8157771759892425, -0.5057053459652427], [0.9731582247379497, 2.2459721722410184, -1.6474522505842983], [-2.011954187959715, 0.6042709716797869, -1.0485616399675888], [-0.40994706586377727, -0.1158729665822918, -2.318109697394429], [1.5770424436171675, 0.0, 0.0], [2.292718146893926, 1.3915527243580657, 0.0], [3.782335574419716, 1.3186147352454611, -0.4807937764724426], [4.419592497958776, 0.3492507324843874, 0.14976318527134924], [3.8473329977166366, 1.1104432852892043, -1.779574085222873], [4.368010411051582, 2.4684340853352, -0.20560554420871258], [2.2938967322202206, 1.8704189044736008, 1.2405689893126395], [1.6494649440008846, 2.2352148943143417, -0.8029305726284784], [1.9974224573334864, -0.6906780683055322, 1.0535722235492928], [1.9277183224308954, -0.6529932317206211, -1.1102241252095348], [-0.35014935725347235, 0.7132786644586319, 1.0706086973199371], [-0.42668432219275443, -1.2491488329668503, 0.15362238828849778], [-3.035262155625421, 1.7756378008307554, -5.067180175425399]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0424', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
