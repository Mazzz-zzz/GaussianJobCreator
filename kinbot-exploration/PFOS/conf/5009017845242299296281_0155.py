import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0155'
logfile = 'conf/5009017845242299296281_0155.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863864, 0.7718203945763884, 1.1635336229088438], [-0.3976197158559589, 2.318865504557597, 1.2514273698287361], [-1.1233882121466392, 3.1329845708321873, 2.379279124290356], [-1.223597108641382, 2.3151670007983443, 3.7137144948532206], [-0.09347781129710513, 1.6642568079316387, 3.927368976831252], [-1.4669716767588634, 3.1280949290491136, 4.737832185168427], [-2.6155546136485306, 1.0577772469323476, 3.6581389404963227], [-2.3999543975651005, 0.10943549990579163, 4.695730134310813], [-3.818868641292418, 1.7838930810170879, 3.4845073822620445], [-2.2605739139485053, 0.3862677937526728, 2.278426511422851], [-0.43435958530525626, 4.245343424397677, 2.6177930333052033], [-2.3533195674001606, 3.4429499950288984, 1.9898948685020783], [0.9186702629217363, 2.411095160260259, 1.4604772899922662], [-0.6863393184328432, 2.8685142811654005, 0.07804489435392688], [-0.2560445759534373, 0.2684346684142533, 2.3183978473060134], [-2.0076024771874503, 0.5917020341966409, 1.0837604470856876], [1.577042443617162, 0.0, 0.0], [2.292718146893916, 1.3915527243580528, 0.0], [1.6005215470082557, 2.4407219045638966, -0.9357086002340271], [1.3760692761371767, 1.914535025410579, -2.125714308102145], [0.46117091336275706, 2.8519467789832507, -0.41882897761345317], [2.404864907693509, 3.478080011182909, -1.068930340330678], [3.540059223330466, 1.2295174124846031, -0.43152105569273796], [2.3088468039522465, 1.896094738758372, 1.2310220414904727], [1.9974224573334778, -0.6906780683055251, 1.0535722235493081], [1.9277183224308962, -0.6529932317206336, -1.1102241252095235], [-0.35014935725347973, -1.2838136616209446, 0.08241309473865535], [-0.42668432219275443, 0.49153352563553854, -1.1586058166012312], [-1.7302689259614916, -0.410574470868132, 2.419072918152931]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0155', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
