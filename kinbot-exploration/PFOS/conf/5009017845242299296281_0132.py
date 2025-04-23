import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0132'
logfile = 'conf/5009017845242299296281_0132.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863844, 0.6217394783082179, -1.2501828803165003], [-0.39761971585595707, -0.0756648590159535, -2.633910119820667], [-1.1233882121466345, 0.49402387891337196, -3.9028837901505375], [-0.3678111583333658, 0.11496132015470119, -5.223868290930238], [-1.1821704538825692, 0.23639236066085972, -6.257587490341371], [0.6879514164795761, 0.9051956426166068, -5.396480983918136], [0.26055307355107626, -1.6528948974501925, -5.185059566666356], [0.563863206493521, -2.046471686712383, -6.517550197784828], [1.1717646888843918, -1.7417603337975769, -4.104710966467124], [-1.075302365678474, -2.3594548912916498, -4.741431512557303], [-2.351535285418468, -0.013330338437763909, -3.9597873669787513], [-1.1885301103712396, 1.817110520581263, -3.82786330478367], [0.9186702629217401, 0.059262854653440736, -2.8183083047232302], [-0.6863393184328435, -1.366668279436511, -2.523228685784665], [-0.25604457595343205, 1.8735740976390451, -1.3916701657561907], [-2.007602477187447, 0.6427130616946762, -1.0543092166280588], [1.5770424436171655, 0.0, 0.0], [2.292718146893913, 1.3915527243580605, 0.0], [1.6005215470082441, 2.4407219045638984, -0.9357086002340271], [1.376069276137157, 1.914535025410574, -2.125714308102145], [0.4611709133627462, 2.85194677898325, -0.4188289776134485], [2.404864907693489, 3.478080011182916, -1.0689303403306853], [3.5400592233304633, 1.2295174124846195, -0.4315210556927483], [2.3088468039522487, 1.89609473875838, 1.2310220414904673], [1.9974224573334856, -0.6906780683055227, 1.053572223549298], [1.927718322430895, -0.6529932317206208, -1.1102241252095322], [-0.35014935725347496, 0.7132786644586343, 1.0706086973199347], [-0.42668432219275365, -1.2491488329668519, 0.1536223882885025], [-0.8827468948257032, -3.1043227962604734, -4.154868147277791]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0132', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
