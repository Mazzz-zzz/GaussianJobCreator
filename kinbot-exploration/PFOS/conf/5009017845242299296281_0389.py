import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0389'
logfile = 'conf/5009017845242299296281_0389.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863808, 0.6217394783082131, -1.2501828803165032], [-2.2709622836291916, 0.6501421835576476, -1.2334320314121792], [-2.9970239643018988, 1.3180406141844374, -2.453401484532624], [-4.47595427893112, 0.8159291580463145, -2.595649047209704], [-4.4929097750280675, -0.3734363274303271, -3.17151725798889], [-5.055401694144313, 0.7379651761413913, -1.4010518970869976], [-5.510548804511903, 1.981353038699881, -3.6411463820487984], [-4.736207008264467, 2.37475269635332, -4.767229194168959], [-6.8059689138757, 1.4158722103092816, -3.727604144144427], [-5.583541308689301, 3.188261715905348, -2.6318960815637973], [-3.0150768858013444, 2.635812688105747, -2.2741184520826354], [-2.351095681931018, 1.0276075788965995, -3.5753886090291105], [-2.6427825570536134, -0.6321679360904253, -1.1859143708981001], [-2.667789377892035, 1.2545909529572552, -0.1201175373359528], [-0.371045161828279, -0.11265730320380458, -2.3156464312139007], [-0.24552532002049007, 1.8598848945507176, -1.4267659957399799], [1.5770424436171666, 0.0, 0.0], [2.2927181468939164, 1.3915527243580583, 0.0], [1.6005215470082508, 2.440721904563894, -0.9357086002340316], [1.3760692761371687, 1.914535025410581, -2.1257143081021432], [0.4611709133627562, 2.85194677898325, -0.41882897761345295], [2.4048649076934985, 3.4780800111829078, -1.0689303403306878], [3.540059223330465, 1.2295174124846102, -0.4315210556927456], [2.308846803952253, 1.8960947387583769, 1.2310220414904716], [1.9974224573334824, -0.690678068305529, 1.0535722235493024], [1.9277183224308967, -0.652993231720628, -1.1102241252095306], [-0.35014935725347734, 0.7132786644586331, 1.070608697319933], [-0.42668432219275454, -1.2491488329668494, 0.15362238828849994], [-5.456557680368095, 2.876135853567083, -1.725025092451242]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0389', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
