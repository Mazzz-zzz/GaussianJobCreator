import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0087'
logfile = 'conf/5009017845242299296281_0087.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863838, 0.771820394576385, 1.163533622908845], [-0.3976197158559582, 2.318865504557592, 1.251427369828744], [-1.1233882121466388, 3.1329845708321784, 2.3792791242903664], [-2.5854106473305416, 3.5221161086622663, 1.966178404989463], [-2.5604417585430426, 4.567045136896721, 1.1574746127251334], [-3.1805547950289808, 2.5023425851706373, 1.3539975127698698], [-3.6360803989189696, 3.9668090929445703, 3.4561463757667554], [-2.8559764250854114, 4.779895070342085, 4.3237637978939345], [-4.91376663469693, 4.3293967871434065, 2.9648167322817356], [-3.7580301869750024, 2.5284120810761657, 4.08597565045167], [-1.180898235093625, 2.3869269963298168, 3.4788618839288983], [-0.4526295328666376, 4.249846884503492, 2.6305382892183653], [0.9186702629217388, 2.411095160260257, 1.4604772899922709], [-0.6863393184328441, 2.8685142811653983, 0.07804489435393573], [-0.256044575953429, 0.2684346684142472, 2.3183978473060147], [-2.0076024771874477, 0.5917020341966329, 1.0837604470856936], [1.577042443617165, 0.0, 0.0], [2.2927181468939155, 1.3915527243580559, 0.0], [3.7823355744197142, 1.3186147352454658, -0.4807937764724475], [4.419592497958773, 0.3492507324843944, 0.14976318527134436], [3.8473329977166335, 1.1104432852891946, -1.7795740852228814], [4.368010411051577, 2.4684340853352102, -0.20560554420872035], [2.2938967322202144, 1.8704189044736115, 1.2405689893126282], [1.649464944000881, 2.2352148943143386, -0.8029305726284901], [1.9974224573334838, -0.6906780683055236, 1.0535722235493006], [1.9277183224308956, -0.6529932317206295, -1.1102241252095286], [-0.35014935725347235, -1.2838136616209457, 0.08241309473865084], [-0.4266843221927558, 0.49153352563553865, -1.1586058166012285], [-4.6307486249364445, 2.409711147404689, 4.486269296521236]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0087', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
