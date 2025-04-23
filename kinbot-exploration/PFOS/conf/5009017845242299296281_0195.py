import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0195'
logfile = 'conf/5009017845242299296281_0195.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863829, 0.7718203945763835, 1.163533622908845], [-0.3976197158559581, 2.318865504557593, 1.2514273698287406], [1.0879414097563231, 2.7598675663867502, 1.4967937649483756], [1.4730420899455323, 2.6857139528931366, 3.015274540831183], [1.0004737759019124, 3.7438354477915525, 3.6505382920974645], [0.9773040172258539, 1.5808858197733988, 3.565273268171844], [3.3329235748106685, 2.64396597813703, 3.2617726137944754], [3.609671411865194, 2.977544939593382, 4.616234395820582], [3.8065427668236014, 1.4737057508841054, 2.6202645210045707], [3.7106516445294835, 3.8726558989253803, 2.3517257942798846], [1.894753649284139, 1.9494330700686116, 0.8176461448894706], [1.2598103705638528, 4.009816458287681, 1.0862835586981863], [-0.7843720924895168, 2.8158789821420873, 0.07319662459686571], [-1.1581983763610575, 2.837834268725449, 2.207717463520164], [-0.25604457595342994, 0.26843466841424707, 2.318397847306014], [-2.0076024771874463, 0.5917020341966348, 1.0837604470856903], [1.5770424436171657, 0.0, 0.0], [2.2927181468939155, 1.3915527243580554, 0.0], [3.782335574419715, 1.3186147352454642, -0.48079377647244664], [4.419592497958774, 0.3492507324843942, 0.14976318527135574], [3.8473329977166375, 1.110443285289196, -1.7795740852228752], [4.368010411051578, 2.468434085335209, -0.20560554420872224], [2.2938967322202144, 1.8704189044736088, 1.2405689893126295], [1.649464944000881, 2.2352148943143355, -0.8029305726284925], [1.9974224573334833, -0.6906780683055245, 1.0535722235493015], [1.927718322430898, -0.65299323172063, -1.110224125209529], [-0.35014935725347274, -1.2838136616209446, 0.08241309473865308], [-0.4266843221927575, 0.4915335256355402, -1.15860581660123], [3.026496535871437, 4.018818008539789, 1.6834961842608678]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0195', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
