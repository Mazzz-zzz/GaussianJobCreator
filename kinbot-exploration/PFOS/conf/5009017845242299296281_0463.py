import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0463'
logfile = 'conf/5009017845242299296281_0463.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863799, 0.6217394783082084, -1.2501828803165038], [-0.34660204151390406, 2.127818130564322, -1.566386336981125], [1.1624535490467058, 2.491486091410015, -1.7945884035232629], [1.9172819807099777, 1.3689994871066282, -2.5882149031527923], [2.227217092188635, 0.37480551409109647, -1.774568472711842], [1.1566135117755332, 0.9104732696570285, -3.5781832000357263], [3.5062602700141103, 2.00649169213748, -3.3565995826981427], [4.313987042727912, 0.8877400152155225, -3.7006971761736698], [3.1421067329679917, 3.0310533505600152, -4.26380539242121], [4.116322039547438, 2.695878460964219, -2.0786664271185944], [1.2315520144708516, 3.6207524606765054, -2.4938554639699166], [1.768541893208141, 2.6511749791394865, -0.6251499999737746], [-1.0117971004258954, 2.410996235613986, -2.6898422467352874], [-0.8215282967123818, 2.883031579863845, -0.5833314461381616], [-2.011954187959715, 0.6042709716797761, -1.0485616399675926], [-0.40994706586377555, -0.11587296658230478, -2.318109697394431], [1.5770424436171662, 0.0, 0.0], [2.292718146893914, 1.3915527243580592, 0.0], [2.341079856722316, 2.059852692894931, 1.416502376706475], [1.1453916410070264, 2.0331725558915923, 1.9759511228307884], [3.2022179922846394, 1.4422259961598312, 2.1984030628363276], [2.721555559149491, 3.3152059779288, 1.274535884539409], [1.629232339193971, 2.2122558673105788, -0.8090479336198866], [3.54558683009438, 1.2600392214310832, -0.4280914688619797], [1.9974224573334867, -0.6906780683055218, 1.0535722235493017], [1.9277183224308967, -0.6529932317206242, -1.1102241252095286], [-0.3501493572534749, 0.7132786644586344, 1.0706086973199327], [-0.426684322192755, -1.2491488329668505, 0.15362238828850233], [4.593016986767815, 3.500540890668066, -2.3261390280538508]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0463', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
