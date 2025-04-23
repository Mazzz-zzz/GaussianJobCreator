import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0299'
logfile = 'conf/5009017845242299296281_0299.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863827, 0.771820394576382, 1.1635336229088498], [-0.34660204151390067, 0.292621294684338, 2.6259377241923803], [-0.7363023803695165, -1.175309686243212, 3.0195004580111453], [-2.245266873553545, -1.285267541889906, 3.4324418428267744], [-2.6409584449074375, -2.543912097586819, 3.356652622940251], [-2.4192802354782867, -0.840213353043175, 4.673629241818011], [-3.350749082422242, -0.2613798489209822, 2.31385952168248], [-4.68981527198154, -0.7105559820462194, 2.4799123723650127], [-2.943049135445392, 1.0873335230058614, 2.4556684257358037], [-2.812813595393545, -0.7787828054582919, 0.926970116513251], [-0.5266683954723089, -1.9706918044823225, 1.9743151008291857], [0.004125934750531798, -1.5779705333152394, 4.044224788693656], [0.9809271994195113, 0.4056495998098974, 2.7245186729056177], [-0.9145558447621658, 1.1282242700023666, 3.4870389756753246], [-2.011954187959716, 0.6059455318059185, 1.047594832227998], [-0.4099470658637743, 2.0654783699937833, 1.0587059160250798], [1.5770424436171664, 0.0, 0.0], [2.2927181468939146, 1.3915527243580534, 0.0], [3.7823355744197107, 1.3186147352454585, -0.480793776472449], [4.419592497958775, 0.34925073248439054, 0.14976318527134785], [3.8473329977166326, 1.1104432852892017, -1.7795740852228814], [4.3680104110515785, 2.468434085335206, -0.20560554420871718], [2.2938967322202197, 1.8704189044736017, 1.2405689893126348], [1.6494649440008835, 2.235214894314338, -0.802930572628484], [1.997422457333482, -0.6906780683055284, 1.053572223549295], [1.927718322430894, -0.6529932317206257, -1.1102241252095342], [-0.35014935725347685, -1.2838136616209446, 0.0824130947386485], [-0.4266843221927575, 0.4915335256355498, -1.158605816601226], [-3.372354856186968, -1.4950326292059495, 0.595515513064823]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0299', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
