import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0483'
logfile = 'conf/5009017845242299296281_0483.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863818, 0.6217394783082144, -1.2501828803165032], [-0.39761971585595246, -0.07566485901595822, -2.6339101198206682], [-0.749183351766654, -1.598819837150671, -2.7680043658220974], [-2.0819874597153207, -1.9555730430184117, -2.022386823820668], [-2.559930149196675, -3.099663561422626, -2.47969373864473], [-1.8676581621436081, -2.056568876678193, -0.71366117100462], [-3.3945978505277346, -0.6377300155556188, -2.2713989681308866], [-3.0839451388502597, 0.46439975006774536, -1.4280509496566538], [-3.6119681777355006, -0.5323922951478793, -3.6667522843762335], [-4.617765080581883, -1.4157570635800159, -1.655705481342929], [-0.8953137364216012, -1.8945036554923982, -4.056495166238402], [0.2259041467951674, -2.3330690032916617, -2.2478963639045073], [-1.1305489784219283, 0.5952337900286001, -3.526959494450473], [0.887865254859783, 0.09037708117406056, -2.9207161247297293], [-0.25604457595342706, 1.8735740976390423, -1.3916701657561992], [-2.0076024771874468, 0.6427130616946742, -1.054309216628068], [1.5770424436171655, 0.0, 0.0], [2.292718146893915, 1.3915527243580577, 0.0], [2.34107985672232, 2.059852692894931, 1.4165023767064724], [1.145391641007031, 2.033172555891597, 1.9759511228307884], [3.202217992284644, 1.4422259961598303, 2.198403062836327], [2.721555559149497, 3.315205977928799, 1.2745358845394046], [1.6292323391939765, 2.2122558673105788, -0.8090479336198869], [3.5455868300943774, 1.2600392214310796, -0.42809146886198113], [1.9974224573334833, -0.6906780683055244, 1.0535722235492997], [1.9277183224308974, -0.6529932317206242, -1.1102241252095306], [-0.35014935725347507, 0.7132786644586373, 1.070608697319929], [-0.4266843221927575, -1.2491488329668499, 0.1536223882885011], [-5.428043559604871, -1.2062048426275314, -2.1409879143879302]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0483', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
