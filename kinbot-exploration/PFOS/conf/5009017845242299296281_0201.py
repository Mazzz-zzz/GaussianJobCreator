import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0201'
logfile = 'conf/5009017845242299296281_0201.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586385, -1.3935598728845984, 0.08664925740765084], [-0.39761971585595673, -2.2432006455416373, 1.3824827499919183], [-1.1233882121466356, -3.6270084497455453, 1.5236046658601603], [-1.2235971086413773, -4.373754595344558, 0.14813618926815572], [-1.448535500070832, -5.66013480745485, 0.350943026009932], [-2.2078139831012678, -3.8623016377645314, -0.5857800571551395], [0.3571901289937569, -4.211693885893079, -0.850117970200745], [0.34893161476605494, -5.211323566125522, -1.8616036632044208], [0.5452803995777125, -2.830743976550133, -1.1010286508131113], [1.3787336598953412, -4.636059172755695, 0.27100001983783], [-0.434359585305253, -4.389746980891082, 2.367678736664982], [-2.3533195674001575, -3.444774504497561, 1.98673472540347], [0.9186702629217397, -2.4703580149136926, 1.35783101473095], [-0.686339318432842, -1.5018460017288895, 2.445183791430722], [-0.25604457595343166, -2.1420087660532854, -0.9267276815498239], [-2.0076024771874486, -1.2344150958913085, -0.02945123045763407], [1.5770424436171633, 0.0, 0.0], [2.2927181468939146, 1.3915527243580548, 0.0], [3.782335574419715, 1.3186147352454565, -0.4807937764724477], [4.419592497958776, 0.34925073248439253, 0.14976318527134747], [3.8473329977166375, 1.1104432852891941, -1.7795740852228803], [4.368010411051579, 2.468434085335206, -0.20560554420871713], [2.29389673222022, 1.8704189044736033, 1.2405689893126275], [1.6494649440008842, 2.2352148943143386, -0.8029305726284843], [1.9974224573334824, -0.6906780683055282, 1.0535722235492944], [1.9277183224308947, -0.6529932317206281, -1.1102241252095348], [-0.3501493572534751, 0.5705349971623106, -1.1530217920585808], [-0.4266843221927543, 0.757615307331305, 1.0049834283127268], [2.1984107742981487, -4.1288736170025055, 0.188207345913032]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0201', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
