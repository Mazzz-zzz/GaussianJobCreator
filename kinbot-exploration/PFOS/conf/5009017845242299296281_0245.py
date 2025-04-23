import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0245'
logfile = 'conf/5009017845242299296281_0245.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863835, 0.6217394783082121, -1.2501828803165012], [-0.39761971585595535, -0.07566485901595937, -2.6339101198206674], [1.0879414097563276, -0.08367235852192167, -3.138512306045848], [1.4730420899455359, 1.268447375297704, -3.833533780919367], [0.9614495666061568, 2.2823561439552384, -3.1576905727862905], [2.7953494776680796, 1.3981654073395322, -3.8933371678685145], [0.8248630198382841, 1.3578016551479501, -5.592380856950237], [0.8535738017776328, 2.7173970689148845, -6.008442162137757], [1.4300365113903912, 0.29621478626472064, -6.308031006524267], [-0.6687885391477252, 0.948893903962471, -5.305427722853217], [1.8947536492841413, -0.26661420225360816, -2.0970816341016434], [1.2598103705638581, -1.0641590715978395, -4.015744696739174], [-0.7843720924895141, -1.34454935469888, -2.4752210448161533], [-1.158198376361053, 0.49302227342428917, -3.5614953002063574], [-0.2560445759534298, 1.8735740976390411, -1.3916701657561943], [-2.007602477187446, 0.6427130616946719, -1.05430921662806], [1.5770424436171646, 0.0, 0.0], [2.292718146893913, 1.3915527243580585, 0.0], [1.6005215470082437, 2.4407219045639037, -0.9357086002340332], [1.3760692761371573, 1.9145350254105766, -2.1257143081021455], [0.46117091336274796, 2.8519467789832476, -0.41882897761345095], [2.404864907693493, 3.478080011182914, -1.0689303403306858], [3.5400592233304655, 1.2295174124846147, -0.431521055692744], [2.308846803952248, 1.896094738758376, 1.2310220414904707], [1.9974224573334847, -0.6906780683055224, 1.0535722235493012], [1.9277183224308976, -0.6529932317206233, -1.1102241252095313], [-0.3501493572534779, 0.7132786644586347, 1.0706086973199314], [-0.4266843221927516, -1.249148832966852, 0.15362238828850366], [-1.2203053177951182, 1.7377083148430834, -5.2076348051470935]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0245', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
