import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0174'
logfile = 'conf/5009017845242299296281_0174.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863817, -1.393559872884596, 0.08664925740765565], [-0.3466020415139008, -2.4204394252486683, -1.0595513872112563], [-0.736302380369517, -2.0273092602547993, -2.5275982746061056], [-0.5090926255515532, -0.4996509844589478, -2.8000337324856055], [-1.5087734085102669, 0.20040764265050626, -2.2930246656340585], [0.6350825859976149, -0.0971565751238951, -2.254358576111915], [-0.41185458969835453, -0.12571598289924055, -4.636438020791137], [-1.4318479285867929, -0.865658559745446, -5.295626597302747], [-0.2399351960211971, 1.2742571457195522, -4.762436960750947], [0.977471611715914, -0.8053695259577733, -4.933736519409601], [0.018370741190138658, -2.722774439312998, -3.373640954864588], [-2.0169579203546326, -2.3003854882603787, -2.7415622040257057], [0.9809271994195118, -2.562327183726283, -1.0109564779824414], [-0.9145558447621661, -3.583976471922484, -0.766448608849453], [-2.0119541879597143, -1.2102165034857024, 0.0009668077396031383], [-0.4099470658637756, -1.9496054034114785, 1.2594037813693577], [1.5770424436171668, 0.0, 0.0], [2.2927181468939137, 1.3915527243580577, 0.0], [3.7823355744197116, 1.3186147352454634, -0.4807937764724417], [4.419592497958777, 0.34925073248440097, 0.1497631852713459], [3.847332997716636, 1.110443285289201, -1.7795740852228776], [4.368010411051578, 2.4684340853352085, -0.2056055442087194], [2.2938967322202157, 1.8704189044736073, 1.2405689893126328], [1.6494649440008775, 2.235214894314339, -0.8029305726284858], [1.9974224573334851, -0.6906780683055247, 1.0535722235492992], [1.9277183224308956, -0.6529932317206254, -1.1102241252095308], [-0.35014935725347646, 0.5705349971623078, -1.1530217920585801], [-0.42668432219275815, 0.7576153073313018, 1.00498342831273], [0.8475108863146382, -1.7056206378326204, -5.263312811210568]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0174', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
