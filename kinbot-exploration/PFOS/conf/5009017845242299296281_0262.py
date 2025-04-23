import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0262'
logfile = 'conf/5009017845242299296281_0262.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863834, 0.6217394783082153, -1.2501828803164994], [-0.3976197158559578, -0.07566485901595062, -2.6339101198206682], [-1.123388212146638, 0.49402387891337435, -3.9028837901505358], [-2.5854106473305416, -0.05829760723787371, -4.033331227674647], [-2.5604417585430426, -1.2811201495928228, -4.53391441514531], [-3.1805547950289808, -0.07857504986565651, -2.8440910041143326], [-3.636080398918971, 1.0097060141392613, -5.163430634336477], [-4.7806825965319915, 0.25837444922166963, -5.547901304983696], [-3.692686028550908, 2.2973105807824634, -4.576721972135597], [-2.647865253360702, 1.0864852305653665, -6.3873907676019535], [-1.180898235093625, 1.819319269574921, -3.8065703577649486], [-0.45262953286663754, 0.15318954183903194, -4.995744508783357], [0.9186702629217391, 0.0592628546534395, -2.818308304723231], [-0.6863393184328441, -1.3666682794365086, -2.523228685784666], [-0.25604457595342744, 1.873574097639046, -1.3916701657561916], [-2.0076024771874477, 0.6427130616946786, -1.0543092166280574], [1.5770424436171655, 0.0, 0.0], [2.292718146893916, 1.3915527243580548, 0.0], [3.7823355744197142, 1.3186147352454634, -0.48079377647244886], [4.419592497958776, 0.3492507324843952, 0.14976318527134774], [3.8473329977166375, 1.1104432852891992, -1.7795740852228792], [4.3680104110515785, 2.468434085335209, -0.20560554420871907], [2.2938967322202215, 1.8704189044736086, 1.2405689893126328], [1.6494649440008806, 2.235214894314339, -0.8029305726284857], [1.9974224573334862, -0.6906780683055294, 1.053572223549296], [1.927718322430895, -0.6529932317206211, -1.1102241252095362], [-0.3501493572534723, 0.7132786644586345, 1.0706086973199351], [-0.426684322192756, -1.2491488329668468, 0.1536223882885019], [-2.701373641589533, 1.95583644087931, -6.808493738632597]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0262', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
