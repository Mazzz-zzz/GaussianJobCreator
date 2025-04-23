import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0482'
logfile = 'conf/5009017845242299296281_0482.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863838, 0.7718203945763837, 1.1635336229088455], [-0.3466020415139057, 0.29262129468434483, 2.625937724192378], [1.1624535490467036, 0.30841610108309786, 3.0549844500983023], [1.322955526797157, 0.3508819546540591, 4.614478254324753], [1.1329590874195934, 1.5833432251349562, 5.0520196956368455], [0.44565046019380916, -0.4645812287616045, 5.192626727607082], [3.0339402189610625, -0.19904667675480195, 5.154638614451783], [3.2390234210468942, 0.2358539450551379, 6.493008270855253], [3.1906412013888583, -1.5386168812243965, 4.72295210661136], [3.881215795346239, 0.7145411838722207, 4.191298618196506], [1.7500514572469126, -0.7950580504033852, 2.6010405683398576], [1.7633282149005816, 1.3769084895275212, 2.5473394429719387], [-1.0117971004259017, 1.123973600038376, 3.4329051118380045], [-0.8215282967123827, -0.9363359387499617, 2.7884443111439525], [-2.011954187959717, 0.6059455318059163, 1.0475948322279898], [-0.40994706586377744, 2.065478369993787, 1.0587059160250722], [1.5770424436171688, 0.0, 0.0], [2.2927181468939133, 1.3915527243580585, 0.0], [3.7823355744197156, 1.3186147352454691, -0.4807937764724413], [4.419592497958776, 0.34925073248439886, 0.1497631852713544], [3.8473329977166397, 1.110443285289201, -1.7795740852228743], [4.368010411051577, 2.468434085335214, -0.2056055442087154], [2.2938967322202135, 1.870418904473612, 1.2405689893126337], [1.649464944000881, 2.2352148943143386, -0.8029305726284872], [1.9974224573334847, -0.6906780683055238, 1.053572223549301], [1.9277183224308987, -0.6529932317206265, -1.1102241252095286], [-0.3501493572534733, -1.2838136616209448, 0.08241309473865079], [-0.42668432219275454, 0.49153352563554814, -1.1586058166012294], [4.1480403748971355, 1.5257812663839516, 4.64589355764243]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0482', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
