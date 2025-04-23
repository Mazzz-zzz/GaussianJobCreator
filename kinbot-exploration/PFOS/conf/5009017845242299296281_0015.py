import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0015'
logfile = 'conf/5009017845242299296281_0015.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863838, 0.6217394783082186, -1.2501828803165005], [-0.3976197158559571, -0.0756648590159506, -2.633910119820668], [-1.1233882121466345, 0.49402387891337546, -3.902883790150537], [-0.3678111583333658, 0.11496132015470467, -5.223868290930238], [0.05889782664840622, -1.1338827048026743, -5.154514224423079], [-1.17364204803218, 0.25007931135430234, -6.273287056826108], [1.1196427974117358, 1.2206221069089234, -5.518246084398576], [1.8011049074626504, 1.3880747859986626, -4.281314396454323], [1.7377767244349351, 0.7696493936980318, -6.709837806361299], [0.3498157742865964, 2.5527005540479863, -5.855108066437026], [-2.3515352854184677, -0.01333033843775951, -3.9597873669787504], [-1.1885301103712391, 1.8171105205812672, -3.82786330478367], [0.91867026292174, 0.05926285465344072, -2.81830830472323], [-0.6863393184328437, -1.3666682794365066, -2.523228685784665], [-0.25604457595342844, 1.8735740976390454, -1.3916701657561923], [-2.0076024771874477, 0.6427130616946822, -1.0543092166280612], [1.5770424436171657, 0.0, 0.0], [2.2927181468939164, 1.3915527243580565, 0.0], [2.341079856722325, 2.0598526928949314, 1.4165023767064704], [1.1453916410070375, 2.033172555891596, 1.9759511228307938], [3.202217992284648, 1.4422259961598296, 2.198403062836322], [2.7215555591495035, 3.3152059779287986, 1.2745358845394015], [1.6292323391939756, 2.2122558673105797, -0.8090479336198926], [3.545586830094382, 1.2600392214310796, -0.4280914688619866], [1.9974224573334816, -0.6906780683055229, 1.0535722235493], [1.9277183224308954, -0.6529932317206202, -1.110224125209534], [-0.35014935725347335, 0.7132786644586356, 1.0706086973199347], [-0.42668432219275776, -1.2491488329668474, 0.15362238828850197], [0.26908555978531656, 3.1034787829301487, -5.063850610388737]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0015', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
