import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0364'
logfile = 'conf/5009017845242299296281_0364.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863853, -1.3935598728845977, 0.08664925740764831], [-0.3976197158559583, -2.2432006455416347, 1.3824827499919197], [-1.123388212146639, -3.6270084497455426, 1.5236046658601627], [-1.22359710864138, -4.373754595344556, 0.14813618926815955], [-2.214848248727895, -3.8691855932816606, -0.5655396772927306], [-0.08476421785157562, -4.253251034978852, -0.5282058686785116], [-1.5460987181777943, -6.208748775927685, 0.37275230664947895], [-0.3251967235864591, -6.8317887185844715, 0.7517649149715777], [-2.756206431121199, -6.3277159809599794, 1.0985917509136176], [-1.853195658006872, -6.571232973367104, -1.1288599289038037], [-0.4343595853052585, -4.389746980891078, 2.3676787366649856], [-2.3533195674001606, -3.4447745044975564, 1.9867347254034717], [0.9186702629217383, -2.470358014913693, 1.3578310147309502], [-0.6863393184328439, -1.5018460017288873, 2.4451837914307233], [-0.2560445759534327, -2.1420087660532854, -0.9267276815498239], [-2.00760247718745, -1.2344150958913065, -0.029451230457633475], [1.5770424436171644, 0.0, 0.0], [2.2927181468939164, 1.391552724358053, 0.0], [2.341079856722326, 2.059852692894924, 1.4165023767064746], [1.145391641007041, 2.0331725558915927, 1.9759511228307973], [3.202217992284648, 1.4422259961598214, 2.19840306283633], [2.7215555591495093, 3.3152059779287937, 1.2745358845394046], [1.6292323391939763, 2.2122558673105748, -0.8090479336198805], [3.545586830094379, 1.2600392214310754, -0.42809146886198185], [1.9974224573334818, -0.6906780683055287, 1.0535722235492955], [1.927718322430895, -0.65299323172063, -1.110224125209535], [-0.3501493572534759, 0.5705349971623116, -1.1530217920585812], [-0.4266843221927576, 0.7576153073313063, 1.0049834283127286], [-1.0558503606080905, -6.9041262491538316, -1.5640437180642348]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0364', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
