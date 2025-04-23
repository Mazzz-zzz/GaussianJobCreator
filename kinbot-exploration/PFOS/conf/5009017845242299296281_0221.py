import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0221'
logfile = 'conf/5009017845242299296281_0221.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863832, 0.6217394783082147, -1.2501828803165012], [-0.39761971585595596, -0.07566485901595645, -2.6339101198206674], [1.087941409756327, -0.0836723585219196, -3.1385123060458486], [1.4730420899455359, 1.2684473752977057, -3.8335337809193666], [2.7879690257787098, 1.3947914162119794, -3.872344601869518], [0.9875700214111515, 1.3076225145166933, -5.071145278802742], [0.7750378848301915, 2.7483816085676356, -2.914704021958527], [1.4740930542051456, 3.909506777415022, -3.3454489815915824], [-0.6348786909392101, 2.6177978171901373, -2.935677682777832], [1.2673983840267085, 2.3703777326979893, -1.4672107271666357], [1.8947536492841406, -0.26661420225360666, -2.097081634101643], [1.2598103705638566, -1.0641590715978364, -4.015744696739173], [-0.7843720924895167, -1.3445493546988778, -2.4752210448161547], [-1.158198376361053, 0.4930222734242932, -3.561495300206358], [-0.25604457595342917, 1.8735740976390423, -1.3916701657561925], [-2.0076024771874468, 0.6427130616946735, -1.0543092166280628], [1.5770424436171653, 0.0, 0.0], [2.2927181468939146, 1.3915527243580563, 0.0], [3.7823355744197147, 1.3186147352454636, -0.4807937764724431], [4.419592497958776, 0.3492507324843961, 0.14976318527134985], [3.8473329977166344, 1.1104432852892008, -1.7795740852228752], [4.368010411051574, 2.468434085335211, -0.2056055442087144], [2.293896732220216, 1.8704189044736084, 1.2405689893126328], [1.6494649440008788, 2.235214894314338, -0.8029305726284861], [1.9974224573334844, -0.6906780683055211, 1.0535722235493001], [1.9277183224308976, -0.6529932317206216, -1.1102241252095313], [-0.3501493572534749, 0.7132786644586353, 1.0706086973199322], [-0.4266843221927544, -1.2491488329668525, 0.15362238828850033], [2.1035634440627096, 2.815156926354197, -1.2698231693988544]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0221', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
