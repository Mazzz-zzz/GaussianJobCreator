import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0039'
logfile = 'conf/5009017845242299296281_0039.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863851, 0.7718203945763887, 1.163533622908844], [-0.3976197158559555, 2.318865504557596, 1.2514273698287408], [-1.1233882121466328, 3.1329845708321837, 2.379279124290362], [-0.36781115833336236, 4.466521985892244, 2.7114935691716706], [0.6688104135818964, 4.21303689453365, 3.491020969700902], [0.05149838568731567, 5.052048948137965, 1.5933532723710033], [-1.4757407605671005, 5.692943503854932, 3.6003386557234682], [-0.6501343745360334, 6.676668829665094, 4.211149253629566], [-2.547144872069914, 5.989869799285394, 2.7232498480279155], [-2.0277986704675164, 4.731543834196867, 4.71903672451002], [-2.351535285418464, 3.435941622607186, 1.9683492717612103], [-1.1885301103712402, 2.4064716038662963, 3.487595524699164], [0.9186702629217429, 2.4110951602602557, 1.4604772899922718], [-0.6863393184328367, 2.868514281165401, 0.07804489435393326], [-0.25604457595343605, 0.2684346684142489, 2.318397847306012], [-2.007602477187449, 0.5917020341966424, 1.0837604470856879], [1.577042443617165, 0.0, 0.0], [2.2927181468939164, 1.3915527243580508, 0.0], [2.3410798567223217, 2.0598526928949297, 1.4165023767064735], [1.1453916410070333, 2.0331725558915945, 1.975951122830791], [3.2022179922846394, 1.4422259961598218, 2.1984030628363302], [2.7215555591495044, 3.3152059779287937, 1.2745358845394041], [1.6292323391939831, 2.2122558673105717, -0.8090479336198882], [3.545586830094382, 1.2600392214310687, -0.4280914688619788], [1.9974224573334782, -0.6906780683055275, 1.0535722235493055], [1.9277183224308982, -0.6529932317206304, -1.1102241252095266], [-0.35014935725347535, -1.2838136616209421, 0.0824130947386495], [-0.4266843221927509, 0.4915335256355413, -1.1586058166012305], [-1.5066768584125763, 4.824100925427159, 5.528868576538308]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0039', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
