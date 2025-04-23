import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0029'
logfile = 'conf/5009017845242299296281_0029.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863841, 0.6217394783082114, -1.250182880316503], [-0.39761971585595607, -0.07566485901596173, -2.6339101198206682], [1.0879414097563267, -0.08367235852192446, -3.138512306045848], [1.3412747433710017, -1.2268379978242032, -4.1818625359304695], [0.31359443784926183, -1.3166825080227558, -5.00787893813258], [2.4483457429192925, -0.9821601217258189, -4.877151381993845], [1.5638064896882031, -2.893937332736271, -3.349442713646216], [2.875485637735201, -2.9464028461138865, -2.8024727050582543], [0.38407981336713964, -3.1435017199249655, -2.6069338871184837], [1.5170785096401858, -3.7993217245233915, -4.637324008175298], [1.3507067731889255, 1.0832029916023942, -3.720189049227472], [1.9101736356829, -0.2683335176786009, -2.113688309236765], [-0.7843720924895156, -1.344549354698883, -2.4752210448161533], [-1.1581983763610542, 0.49302227342428684, -3.5614953002063583], [-0.25604457595343083, 1.8735740976390398, -1.3916701657561963], [-2.0076024771874477, 0.6427130616946711, -1.0543092166280652], [1.5770424436171644, 0.0, 0.0], [2.292718146893912, 1.391552724358055, 0.0], [1.6005215470082452, 2.4407219045638975, -0.9357086002340325], [1.3760692761371616, 1.9145350254105735, -2.125714308102145], [0.4611709133627464, 2.8519467789832476, -0.4188289776134524], [2.404864907693491, 3.4780800111829118, -1.0689303403306898], [3.5400592233304655, 1.22951741248461, -0.4315210556927439], [2.3088468039522465, 1.8960947387583758, 1.2310220414904702], [1.9974224573334831, -0.6906780683055249, 1.0535722235493004], [1.9277183224308956, -0.6529932317206268, -1.1102241252095304], [-0.350149357253478, 0.7132786644586327, 1.0706086973199302], [-0.4266843221927555, -1.2491488329668508, 0.1536223882885046], [1.0330383776779406, -3.3517466648093017, -5.345380576019921]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0029', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
