import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0229'
logfile = 'conf/5009017845242299296281_0229.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863858, -1.3935598728845968, 0.0866492574076495], [-0.39761971585595945, -2.2432006455416365, 1.382482749991915], [-1.1233882121466392, -3.627008449745546, 1.5236046658601567], [-1.223597108641381, -4.373754595344557, 0.1481361892681518], [-0.09347781129710421, -4.23332970793658, -0.5223958143256578], [-1.4669716767588619, -5.667130495747955, 0.3400935814215636], [-2.615554613648528, -3.6969298765090755, -0.9130075028596097], [-3.8446139163399393, -4.210509405394753, -0.4149308965155085], [-2.3519919332735246, -2.3191721128894187, -1.107108069441228], [-2.2627755448334397, -4.443427792507679, -2.254212341735725], [-0.43435958530525737, -4.389746980891081, 2.3676787366649816], [-2.353319567400161, -3.444774504497559, 1.986734725403469], [0.9186702629217371, -2.4703580149136943, 1.3578310147309451], [-0.6863393184328461, -1.501846001728892, 2.4451837914307206], [-0.25604457595343627, -2.142008766053284, -0.9267276815498254], [-2.00760247718745, -1.2344150958913065, -0.029451230457633475], [1.5770424436171628, 0.0, 0.0], [2.292718146893914, 1.3915527243580534, 0.0], [3.782335574419718, 1.318614735245452, -0.4807937764724451], [4.419592497958775, 0.3492507324843881, 0.14976318527134785], [3.8473329977166397, 1.1104432852891923, -1.7795740852228819], [4.3680104110515785, 2.468434085335205, -0.20560554420871707], [2.29389673222022, 1.8704189044736013, 1.2405689893126275], [1.6494649440008797, 2.2352148943143404, -0.8029305726284819], [1.9974224573334798, -0.6906780683055302, 1.0535722235492944], [1.9277183224308918, -0.6529932317206284, -1.1102241252095362], [-0.3501493572534786, 0.5705349971623128, -1.1530217920585801], [-0.4266843221927536, 0.7576153073313053, 1.0049834283127275], [-2.4573667611020404, -3.878006383784475, -3.0147353695294803]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0229', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
