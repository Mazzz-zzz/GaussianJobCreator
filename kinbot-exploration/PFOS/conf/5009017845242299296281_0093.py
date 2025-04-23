import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0093'
logfile = 'conf/5009017845242299296281_0093.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863827, 0.771820394576387, 1.1635336229088462], [-0.397619715855956, 2.318865504557596, 1.2514273698287393], [-1.123388212146635, 3.132984570832184, 2.37927912429036], [-0.3678111583333658, 4.466521985892245, 2.7114935691716657], [0.058897826648405066, 5.030881614919977, 1.595285884940587], [-1.1736420480321812, 5.307786300766389, 3.353218565006771], [1.1196427974117347, 4.168630239968733, 3.816212795203277], [0.6538340252901031, 3.9677882071751815, 5.144789354340119], [1.9686989699044857, 3.2677645288266604, 3.1285345578797092], [1.7659468439865353, 5.601582278619653, 3.718836390893924], [-2.3515352854184677, 3.4359416226071815, 1.9683492717612139], [-1.18853011037124, 2.4064716038662977, 3.487595524699166], [0.9186702629217399, 2.411095160260259, 1.4604772899922664], [-0.6863393184328432, 2.868514281165398, 0.0780448943539319], [-0.256044575953429, 0.2684346684142534, 2.3183978473060143], [-2.0076024771874468, 0.5917020341966366, 1.0837604470856936], [1.577042443617166, 0.0, 0.0], [2.292718146893916, 1.3915527243580537, 0.0], [2.3410798567223328, 2.059852692894925, 1.416502376706475], [1.1453916410070426, 2.0331725558915923, 1.9759511228307927], [3.202217992284644, 1.4422259961598194, 2.1984030628363307], [2.721555559149511, 3.315205977928793, 1.2745358845394053], [1.629232339193984, 2.212255867310575, -0.8090479336198857], [3.545586830094387, 1.2600392214310685, -0.42809146886197963], [1.9974224573334842, -0.6906780683055231, 1.0535722235493028], [1.9277183224308974, -0.6529932317206303, -1.110224125209527], [-0.3501493572534729, -1.2838136616209441, 0.08241309473865532], [-0.42668432219275504, 0.49153352563553615, -1.1586058166012283], [1.489920817843396, 6.146429166136023, 4.469112265990112]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0093', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
