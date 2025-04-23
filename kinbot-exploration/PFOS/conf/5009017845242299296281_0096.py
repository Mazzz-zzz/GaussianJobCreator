import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0096'
logfile = 'conf/5009017845242299296281_0096.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863829, -1.3935598728845986, 0.08664925740765458], [-0.34660204151390256, -2.4204394252486705, -1.059551387211251], [-0.9873929842445245, -3.8498102687678504, -0.9708648395635537], [-1.0274356123761041, -4.556709010365347, -2.3702489598329652], [-2.040447113782524, -4.09095957372802, -3.079746742329373], [0.10871386352754484, -4.345207077208824, -3.0285280353229274], [-1.2388794809575765, -6.414242029259144, -2.2073812058925713], [-1.6618516099518503, -6.925713629081326, -3.4651261675121288], [-0.11466179393237848, -6.894962207953295, -1.4929587303054754], [-2.467331563379808, -6.417814460290127, -1.2217670071626487], [-0.2614927667443652, -4.595762031236105, -0.1428452777114611], [-2.230877915023899, -3.758011046649702, -0.5172703827718971], [-0.778453490345139, -1.8458422644301349, -2.185681892820454], [0.9731582247379446, -2.5497215866483653, -1.1213428320614776], [-2.0119541879597187, -1.2102165034856989, 0.0009668077396009857], [-0.40994706586377944, -1.9496054034114783, 1.2594037813693602], [1.577042443617165, 0.0, 0.0], [2.2927181468939186, 1.391552724358055, 0.0], [3.7823355744197187, 1.3186147352454558, -0.480793776472442], [4.419592497958776, 0.34925073248439353, 0.14976318527135296], [3.8473329977166393, 1.110443285289191, -1.7795740852228765], [4.36801041105158, 2.4684340853352023, -0.2056055442087204], [2.2938967322202126, 1.8704189044736101, 1.2405689893126273], [1.6494649440008855, 2.2352148943143337, -0.8029305726284925], [1.9974224573334813, -0.6906780683055226, 1.0535722235493015], [1.9277183224308974, -0.6529932317206297, -1.1102241252095277], [-0.35014935725347585, 0.5705349971623044, -1.1530217920585835], [-0.42668432219275904, 0.7576153073313052, 1.0049834283127235], [-2.3988383294137656, -7.1575762955859545, -0.6020707181334033]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0096', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
