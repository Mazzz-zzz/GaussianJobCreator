import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0049'
logfile = 'conf/5009017845242299296281_0049.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863823, -1.3935598728845977, 0.08664925740765328], [-2.270962283629192, -1.3932545648232157, 0.05367636867327985], [-2.997023964301902, -2.7837283183799233, 0.08524408716294148], [-4.4759542789311215, -2.655862593215633, 0.5912091450482932], [-5.168854707529478, -3.7205321371650713, 0.2266961805165711], [-4.503853254513123, -2.541450369721436, 1.915985175488761], [-5.326480128927349, -1.1401075478729001, -0.11643183180454744], [-6.72993834376042, -1.3000142532958938, 0.04838895750413604], [-4.598236392479623, -0.012179574876579263, 0.3338972919539449], [-4.967238916287429, -1.3650420974921038, -1.6333012491696688], [-3.0150768858013453, -3.2873506947713778, -1.145621521475612], [-2.351095681931021, -3.6101811532690156, 0.8977600360686738], [-2.6427825570536148, -0.7109480038655813, 1.1404306775613342], [-2.6677893778920336, -0.7313203152515864, -1.0264488679511385], [-0.3710451618282827, -1.9490799840121111, 1.25538730210329], [-0.24552532002048985, -2.1655580448419784, -0.8973245689258752], [1.5770424436171662, 0.0, 0.0], [2.292718146893919, 1.3915527243580548, 0.0], [3.782335574419717, 1.318614735245459, -0.48079377647244625], [4.419592497958775, 0.3492507324843912, 0.1497631852713449], [3.8473329977166353, 1.1104432852891994, -1.7795740852228805], [4.368010411051584, 2.4684340853352014, -0.2056055442087228], [2.2938967322202215, 1.8704189044736017, 1.2405689893126308], [1.6494649440008824, 2.235214894314335, -0.8029305726284871], [1.9974224573334842, -0.6906780683055261, 1.0535722235493], [1.9277183224308931, -0.6529932317206278, -1.1102241252095348], [-0.3501493572534748, 0.5705349971623102, -1.1530217920585797], [-0.4266843221927533, 0.7576153073313034, 1.0049834283127268], [-4.808842825024435, -0.5166840229162979, -2.070521608108543]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0049', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
