import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0264'
logfile = 'conf/5009017845242299296281_0264.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863873, -1.3935598728845968, 0.08664925740765322], [-0.3466020415139064, -2.420439425248666, -1.0595513872112565], [-0.7363023803695208, -2.027309260254797, -2.5275982746061074], [0.11232051632789275, -2.815384817700021, -3.585124280886944], [1.3083470166553104, -2.265328906921625, -3.701146040093135], [0.2433393664383972, -4.087232389653753, -3.218875929103458], [-0.6978983773361233, -2.8041965493719943, -5.277780033872964], [-1.765667469854746, -3.7434309186721997, -5.26599137842827], [-0.8391756380146922, -1.4471552235726262, -5.657037660475062], [0.4974230088710358, -3.4088054948556596, -6.106150383732213], [-2.0209682759934653, -2.311513884043151, -2.7221562889464304], [-0.5312999584009176, -0.7299509739637096, -2.7153888107570547], [0.980927199419505, -2.562327183726282, -1.0109564779824385], [-0.9145558447621747, -3.583976471922481, -0.7664486088494523], [-2.01195418795972, -1.2102165034856947, 0.0009668077395966828], [-0.4099470658637822, -1.949605403411482, 1.259403781369355], [1.5770424436171635, 0.0, 0.0], [2.292718146893915, 1.3915527243580528, 0.0], [2.3410798567223208, 2.0598526928949275, 1.4165023767064746], [1.1453916410070426, 2.0331725558915954, 1.975951122830793], [3.2022179922846457, 1.4422259961598258, 2.1984030628363254], [2.7215555591495058, 3.315205977928791, 1.274535884539405], [1.6292323391939798, 2.212255867310575, -0.8090479336198887], [3.545586830094382, 1.2600392214310712, -0.42809146886198157], [1.9974224573334811, -0.6906780683055278, 1.0535722235493008], [1.927718322430892, -0.6529932317206268, -1.1102241252095317], [-0.35014935725347623, 0.5705349971623095, -1.1530217920585835], [-0.42668432219276103, 0.7576153073313076, 1.0049834283127246], [0.5355286641010409, -3.0053071733428878, -6.984614640416492]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0264', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
