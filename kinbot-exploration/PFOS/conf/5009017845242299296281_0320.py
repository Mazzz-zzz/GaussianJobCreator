import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0320'
logfile = 'conf/5009017845242299296281_0320.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863844, 0.7718203945763873, 1.1635336229088467], [-0.39761971585595707, 2.3188655045575963, 1.2514273698287368], [-1.1233882121466345, 3.1329845708321864, 2.3792791242903597], [-0.36781115833336575, 4.466521985892248, 2.7114935691716626], [-1.1821704538825695, 5.301033552708925, 3.333515534763536], [0.6879514164795761, 4.220891801804463, 3.4821629138600065], [0.26055307355107626, 5.3168407535936995, 1.1610808123556096], [0.563863206493521, 6.667599885078152, 1.4864786300738655], [1.1717646888843916, 4.425664139051895, 0.5439467868607738], [-1.0753023656784737, 5.285927585824527, 0.3273678813366047], [-2.351535285418468, 3.4359416226071846, 1.9683492717612097], [-1.1885301103712396, 2.4064716038663, 3.4875955246991626], [0.91867026292174, 2.4110951602602593, 1.4604772899922667], [-0.6863393184328435, 2.8685142811653987, 0.07804489435392938], [-0.25604457595343205, 0.26843466841425323, 2.318397847306013], [-2.0076024771874477, 0.5917020341966378, 1.083760447085691], [1.5770424436171655, 0.0, 0.0], [2.2927181468939155, 1.391552724358054, 0.0], [3.7823355744197142, 1.3186147352454565, -0.4807937764724459], [4.419592497958775, 0.34925073248439287, 0.14976318527135224], [3.847332997716637, 1.110443285289188, -1.77957408522288], [4.368010411051578, 2.468434085335205, -0.20560554420872612], [2.2938967322202157, 1.8704189044736086, 1.2405689893126273], [1.6494649440008828, 2.235214894314331, -0.8029305726284939], [1.9974224573334838, -0.6906780683055235, 1.053572223549306], [1.9277183224308976, -0.6529932317206303, -1.1102241252095257], [-0.3501493572534726, -1.2838136616209455, 0.0824130947386554], [-0.42668432219275615, 0.49153352563553454, -1.158605816601228], [-1.6453231695648611, 4.566270298864316, 0.6325433376855504]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0320', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
