import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0330'
logfile = 'conf/5009017845242299296281_0330.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863852, 0.7718203945763825, 1.1635336229088462], [-2.2709622836291934, 0.7431123812655658, 1.1797556627388968], [-2.9970239643019028, 1.4656877041954837, 2.3681573973696772], [-2.2704853996782592, 2.79602583912399, 2.7705593355751965], [-1.213602341516464, 2.522841108563524, 3.5154554675534557], [-1.8869755643483324, 3.4609727513360196, 1.684435102028263], [-3.3935951456019624, 3.934431604618534, 3.7525573716297385], [-4.105849023147379, 3.1480681011553213, 4.69953037132539], [-2.628917663197577, 5.078609303675161, 4.0864473423284355], [-4.377248131761582, 4.34053681237623, 2.5915026817635187], [-4.239914773111564, 1.7601618341371414, 1.9973488528651422], [-3.0230545829338547, 0.6702887991292583, 3.4298036615283594], [-2.6427825570536156, 1.3431159399560069, 0.045483693336762755], [-2.667789377892034, -0.5232706377056721, 1.1465664052870872], [-0.371045161828282, 2.061737287215913, 1.060259129110611], [-0.24552532002049043, 0.3056731502912651, 2.324090564665852], [1.5770424436171642, 0.0, 0.0], [2.2927181468939177, 1.3915527243580579, 0.0], [2.3410798567223225, 2.059852692894928, 1.4165023767064733], [1.1453916410070375, 2.0331725558915945, 1.9759511228307942], [3.202217992284641, 1.4422259961598296, 2.1984030628363285], [2.7215555591495035, 3.315205977928792, 1.2745358845394048], [1.6292323391939787, 2.212255867310575, -0.8090479336198882], [3.5455868300943814, 1.2600392214310703, -0.42809146886198046], [1.9974224573334798, -0.6906780683055257, 1.053572223549294], [1.9277183224308934, -0.6529932317206277, -1.110224125209533], [-0.35014935725347457, -1.283813661620942, 0.08241309473864832], [-0.42668432219275576, 0.49153352563554714, -1.1586058166012299], [-4.655145616327261, 5.26161463684365, 2.693242624963989]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0330', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
