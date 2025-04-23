import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0130'
logfile = 'conf/5009017845242299296281_0130.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863832, 0.6217394783082077, -1.2501828803165065], [-0.34660204151390306, 2.1278181305643176, -1.566386336981132], [-0.7363023803695211, 3.2026189464979993, -0.49190218340505165], [-2.2452668735535517, 3.6152156038456056, -0.6031465794771679], [-2.403843994986881, 4.474711320489256, -1.5944196764116785], [-3.003380405095209, 2.5441701834335055, -0.8200004963445041], [-2.8645095706331465, 4.425819047916744, 0.9720117838909693], [-3.096179841055352, 3.4059529780608777, 1.9356678912079093], [-2.0484924953511667, 5.560777482547193, 1.1986000804939887], [-4.253583546146212, 4.937084584635503, 0.4338448290107826], [-0.5266683954723127, 2.6951529346344714, 0.7195116152968796], [0.004125934750526012, 4.2913866722810745, -0.6555498260725754], [0.9809271994195083, 2.1566775839163737, -1.7135621949231843], [-0.9145558447621682, 2.4557522019200992, -2.7205903668258786], [-2.011954187959715, 0.604270971679775, -1.0485616399676003], [-0.40994706586377383, -0.11587296658231108, -2.3181096973944326], [1.5770424436171653, 0.0, 0.0], [2.2927181468939115, 1.3915527243580557, 0.0], [1.6005215470082383, 2.440721904563894, -0.9357086002340349], [1.3760692761371605, 1.9145350254105669, -2.125714308102148], [0.4611709133627413, 2.8519467789832467, -0.41882897761345517], [2.404864907693484, 3.478080011182908, -1.068930340330692], [3.540059223330463, 1.2295174124846124, -0.4315210556927441], [2.30884680395224, 1.8960947387583809, 1.2310220414904693], [1.9974224573334853, -0.6906780683055201, 1.0535722235493044], [1.9277183224308951, -0.6529932317206308, -1.1102241252095248], [-0.3501493572534805, 0.7132786644586357, 1.0706086973199285], [-0.4266843221927554, -1.2491488329668508, 0.15362238828850236], [-4.476685266384857, 5.788624283307973, 0.835177796392429]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0130', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
