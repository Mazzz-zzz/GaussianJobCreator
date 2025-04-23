import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0419'
logfile = 'conf/5009017845242299296281_0419.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863829, -1.3935598728845988, 0.08664925740765335], [-0.34660204151390317, -2.4204394252486683, -1.0595513872112563], [1.1624535490467067, -2.7999021924931125, -1.2603960465750392], [1.9243606525366, -1.7238039150214746, -2.1095693044418984], [3.2268445905155723, -1.8415003694175174, -1.9193455440795668], [1.6500667572714027, -1.8748031547393536, -3.40218762752798], [1.4242023469536254, 0.021444945715179648, -1.6345836051110882], [2.4092408723075187, 0.9177386934024212, -2.1333697609287694], [0.03808133731720843, 0.1384950436761862, -1.8997714495909097], [1.6020493994449865, -0.10743030067361281, -0.07499897239718098], [1.7462981079693125, -2.886492824359217, -0.0684987582168957], [1.2535978601361857, -3.964626716020742, -1.8893151990408876], [-1.0117971004258968, -3.534969835652362, -0.743062865102721], [-0.8215282967123797, -1.9466956411138838, -2.205112865005794], [-2.0119541879597174, -1.2102165034856966, 0.0009668077395988341], [-0.40994706586377827, -1.949605403411481, 1.2594037813693544], [1.5770424436171664, 0.0, 0.0], [2.292718146893919, 1.3915527243580545, 0.0], [1.6005215470082572, 2.440721904563898, -0.9357086002340306], [1.3760692761371696, 1.9145350254105773, -2.1257143081021455], [0.46117091336276106, 2.8519467789832538, -0.4188289776134483], [2.4048649076935056, 3.478080011182906, -1.0689303403306887], [3.54005922333047, 1.2295174124846004, -0.43152105569274385], [2.308846803952253, 1.8960947387583682, 1.2310220414904751], [1.9974224573334824, -0.6906780683055274, 1.0535722235493017], [1.9277183224308962, -0.6529932317206276, -1.110224125209529], [-0.3501493572534736, 0.5705349971623083, -1.15302179205858], [-0.4266843221927568, 0.7576153073313056, 1.0049834283127257], [0.9264477630970549, 0.41085426669705977, 0.38424597627987317]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0419', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
