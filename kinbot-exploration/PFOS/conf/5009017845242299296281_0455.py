import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0455'
logfile = 'conf/5009017845242299296281_0455.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863818, -1.3935598728845981, 0.08664925740764959], [-2.2709622836291916, -1.3932545648232184, 0.053676368673268804], [-2.9970239643018983, -2.783728318379928, 0.08524408716293172], [-3.0516123514917317, -3.443374587786915, -1.336540729974832], [-1.9012352660876388, -3.260558258158618, -1.9608104585981438], [-3.297380408862279, -4.746387544852177, -1.2333160417972227], [-4.411157127461295, -2.702806947322468, -2.397135328841796], [-4.1521016093511784, -3.0427044060810884, -3.753524271540801], [-5.635920736913452, -2.960615050138283, -1.7346809117440998], [-4.063991433386536, -1.1816382839682047, -2.182508889893573], [-2.3329497198238167, -3.59705751812083, 0.9016340449500703], [-4.241488402883601, -2.631664135309427, 0.5195418500188189], [-2.6427825570536188, -0.7109480038655849, 1.140430677561324], [-2.6677893778920305, -0.7313203152515916, -1.0264488679511485], [-0.3710451618282831, -1.949079984012112, 1.255387302103288], [-0.24552532002048527, -2.1655580448419793, -0.8973245689258766], [1.5770424436171664, 0.0, 0.0], [2.2927181468939146, 1.3915527243580594, 0.0], [3.7823355744197125, 1.3186147352454674, -0.48079377647243837], [4.419592497958776, 0.3492507324844012, 0.1497631852713619], [3.8473329977166415, 1.110443285289199, -1.779574085222871], [4.368010411051581, 2.4684340853352067, -0.20560554420870847], [2.2938967322202117, 1.8704189044736081, 1.2405689893126355], [1.6494649440008828, 2.235214894314336, -0.8029305726284883], [1.9974224573334816, -0.6906780683055227, 1.053572223549305], [1.927718322430898, -0.6529932317206228, -1.110224125209533], [-0.3501493572534704, 0.5705349971623063, -1.1530217920585824], [-0.4266843221927614, 0.7576153073313028, 1.004983428312726], [-4.873758731702885, -0.653108977587345, -2.1525742682989986]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0455', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
