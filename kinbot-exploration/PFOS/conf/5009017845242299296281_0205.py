import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0205'
logfile = 'conf/5009017845242299296281_0205.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863809, 0.7718203945763825, 1.1635336229088504], [-0.3466020415139007, 0.29262129468433806, 2.6259377241923807], [-0.7363023803695176, -1.1753096862432113, 3.0195004580111435], [-2.245266873553548, -1.285267541889908, 3.4324418428267705], [-2.4038439949868793, -0.8565477161783808, 4.672423516351342], [-3.0033804050952075, -0.5619438307665817, 2.61331625857658], [-2.8645095706331407, -3.0546964215858043, 3.3468658361034587], [-4.070659490453873, -3.1404977586803073, 4.095378701242879], [-2.735782524459431, -3.4706499972509, 1.9993612008153936], [-1.7096109057685829, -3.7399422308852626, 4.16978131447326], [-0.5266683954723095, -1.9706918044823236, 1.9743151008291802], [0.004125934750529869, -1.5779705333152427, 4.044224788693654], [0.9809271994195111, 0.40564959980989734, 2.7245186729056172], [-0.9145558447621667, 1.1282242700023666, 3.4870389756753246], [-2.0119541879597156, 0.6059455318059191, 1.0475948322279969], [-0.40994706586377383, 2.065478369993785, 1.058705916025083], [1.5770424436171668, 0.0, 0.0], [2.2927181468939173, 1.3915527243580539, 0.0], [3.7823355744197116, 1.318614735245462, -0.48079377647244487], [4.419592497958776, 0.34925073248439376, 0.14976318527135207], [3.8473329977166357, 1.110443285289205, -1.7795740852228774], [4.36801041105158, 2.4684340853352085, -0.20560554420871652], [2.2938967322202197, 1.8704189044736013, 1.2405689893126368], [1.649464944000882, 2.235214894314338, -0.8029305726284842], [1.997422457333483, -0.690678068305529, 1.0535722235492946], [1.9277183224308942, -0.6529932317206227, -1.1102241252095344], [-0.35014935725347723, -1.2838136616209426, 0.08241309473864608], [-0.42668432219275654, 0.49153352563555214, -1.1586058166012272], [-0.905213310140795, -3.203884028512854, 4.130430348621151]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0205', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
