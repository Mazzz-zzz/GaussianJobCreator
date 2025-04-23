import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0031'
logfile = 'conf/5009017845242299296281_0031.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863808, 0.7718203945763855, 1.1635336229088482], [-0.3466020415139017, 0.2926212946843449, 2.625937724192379], [1.1624535490467078, 0.3084161010830951, 3.0549844500983014], [1.9172819807099792, 1.5569601130304838, 2.479695785178603], [3.0265100569963668, 1.7669301563805528, 3.166738432144982], [2.2192278579630744, 1.3657546484456962, 1.1986145238701023], [0.8705838155292759, 3.111389297067898, 2.578873373045672], [0.23042687358553068, 3.1362809130679214, 3.8485231298463938], [1.6582993017555376, 4.170877779901709, 2.066676888139888], [-0.19649889051412675, 2.737446181345754, 1.482495849344955], [1.2315520144708536, 0.3493659548263219, 4.382591343745828], [1.7685418932081418, -0.7841917084166135, 2.6085598817993603], [-1.0117971004258943, 1.123973600038377, 3.4329051118380085], [-0.8215282967123801, -0.9363359387499605, 2.788444311143957], [-2.011954187959715, 0.6059455318059201, 1.047594832227992], [-0.4099470658637734, 2.0654783699937873, 1.0587059160250738], [1.577042443617167, 0.0, 0.0], [2.2927181468939186, 1.3915527243580543, 0.0], [2.341079856722327, 2.059852692894927, 1.4165023767064757], [1.1453916410070382, 2.033172555891592, 1.9759511228307924], [3.2022179922846483, 1.442225996159825, 2.1984030628363276], [2.7215555591495066, 3.3152059779287932, 1.2745358845394033], [1.6292323391939805, 2.212255867310575, -0.8090479336198875], [3.545586830094382, 1.2600392214310732, -0.4280914688619811], [1.9974224573334836, -0.6906780683055282, 1.0535722235492997], [1.9277183224308967, -0.6529932317206292, -1.1102241252095313], [-0.35014935725347573, -1.2838136616209437, 0.08241309473865072], [-0.42668432219275343, 0.49153352563554625, -1.1586058166012276], [0.1733492503863006, 2.097852260549216, 0.8579205089346824]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0031', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
