import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0188'
logfile = 'conf/5009017845242299296281_0188.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863806, -1.3935598728846001, 0.08664925740765343], [-0.3466020415139019, -2.4204394252486656, -1.0595513872112614], [-0.9873929842445216, -3.849810268767846, -0.9708648395635691], [-2.4621724959026183, -3.793768812064694, -0.4403274647009679], [-3.0952492719626608, -2.767547307945039, -0.9812388836715901], [-3.1051333564033308, -4.918796415434231, -0.7400177449467293], [-2.525270366870939, -3.5965120887216333, 1.4248191179809988], [-3.8302771896164236, -3.1584243629303868, 1.7818448102060611], [-1.8906703727799765, -4.734092627903814, 1.980292012028714], [-1.5307227375406445, -2.382667245020261, 1.5589770034074195], [-0.9969507321948189, -4.392760927004433, -2.184992009459939], [-0.27888546472280534, -4.6100562569462165, -0.14595700608829104], [-0.7784534903451409, -1.8458422644301253, -2.1856818928204627], [0.9731582247379456, -2.5497215866483582, -1.1213428320614907], [-2.011954187959716, -1.2102165034857006, 0.0009668077396009871], [-0.40994706586377194, -1.9496054034114845, 1.2594037813693495], [1.5770424436171664, 0.0, 0.0], [2.292718146893918, 1.3915527243580577, 0.0], [3.7823355744197116, 1.3186147352454674, -0.4807937764724472], [4.419592497958776, 0.349250732484398, 0.14976318527134597], [3.847332997716634, 1.110443285289206, -1.7795740852228827], [4.368010411051576, 2.4684340853352063, -0.2056055442087189], [2.293896732220218, 1.8704189044736046, 1.2405689893126322], [1.649464944000878, 2.2352148943143395, -0.8029305726284882], [1.997422457333487, -0.6906780683055261, 1.0535722235492955], [1.9277183224308925, -0.6529932317206212, -1.1102241252095355], [-0.35014935725347757, 0.5705349971623103, -1.1530217920585812], [-0.4266843221927559, 0.7576153073313018, 1.0049834283127281], [-2.0204583922699877, -1.5489613807935236, 1.5914717281603088]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0188', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
