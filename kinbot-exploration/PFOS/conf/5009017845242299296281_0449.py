import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0449'
logfile = 'conf/5009017845242299296281_0449.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863851, 0.7718203945763864, 1.1635336229088453], [-0.3976197158559567, 2.318865504557595, 1.2514273698287417], [-1.1233882121466328, 3.132984570832182, 2.379279124290365], [-0.3678111583333635, 4.466521985892242, 2.7114935691716724], [0.6688104135818954, 4.213036894533651, 3.4910209697009034], [0.05149838568730977, 5.052048948137963, 1.5933532723710029], [-1.475740760567105, 5.692943503854929, 3.6003386557234704], [-2.3402363179501737, 6.288045526114909, 2.640686623903079], [-1.9223669601590467, 5.0596208829030624, 4.785574609986048], [-0.3736644517404418, 6.742184289182709, 4.006619609912503], [-2.351535285418464, 3.4359416226071833, 1.968349271761215], [-1.1885301103712402, 2.4064716038662963, 3.487595524699164], [0.9186702629217423, 2.411095160260256, 1.460477289992272], [-0.6863393184328396, 2.868514281165399, 0.07804489435393576], [-0.2560445759534336, 0.268434668414247, 2.3183978473060134], [-2.007602477187449, 0.591702034196639, 1.083760447085691], [1.5770424436171657, 0.0, 0.0], [2.2927181468939155, 1.3915527243580517, 0.0], [1.6005215470082543, 2.4407219045638975, -0.9357086002340333], [1.3760692761371667, 1.9145350254105704, -2.1257143081021477], [0.4611709133627524, 2.8519467789832467, -0.4188289776134555], [2.4048649076934994, 3.478080011182908, -1.068930340330688], [3.5400592233304633, 1.2295174124846047, -0.43152105569274124], [2.308846803952249, 1.8960947387583724, 1.2310220414904685], [1.9974224573334807, -0.6906780683055271, 1.0535722235493048], [1.927718322430898, -0.6529932317206291, -1.1102241252095273], [-0.3501493572534731, -1.283813661620944, 0.08241309473865303], [-0.42668432219275243, 0.4915335256355401, -1.158605816601229], [-0.33463778567007546, 7.45743282229004, 3.3563436362526695]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0449', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
