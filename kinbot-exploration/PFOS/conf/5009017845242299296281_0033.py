import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0033'
logfile = 'conf/5009017845242299296281_0033.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863816, -1.393559872884599, 0.08664925740765336], [-0.3466020415139052, -2.420439425248664, -1.0595513872112607], [-0.9873929842445249, -3.8498102687678446, -0.970864839563567], [-2.462172495902622, -3.793768812064692, -0.44032746470096756], [-3.095249271962663, -2.767547307945036, -0.9812388836715891], [-3.1051333564033357, -4.91879641543423, -0.7400177449467268], [-2.5252703668709424, -3.5965120887216298, 1.424819117981001], [-2.239776335865652, -4.8582981765830136, 2.0154290015334233], [-1.8335865604435533, -2.400495679065634, 1.7354234183499955], [-4.066604391334009, -3.3076173168432943, 1.5710103134164815], [-0.9969507321948233, -4.392760927004433, -2.184992009459939], [-0.2788854647228104, -4.610056256946216, -0.145957006088289], [-0.7784534903451428, -1.845842264430125, -2.1856818928204618], [0.9731582247379439, -2.549721586648361, -1.121342832061489], [-2.0119541879597174, -1.210216503485699, 0.0009668077396009859], [-0.40994706586377433, -1.949605403411482, 1.2594037813693504], [1.5770424436171646, 0.0, 0.0], [2.29271814689392, 1.3915527243580554, 0.0], [1.600521547008262, 2.4407219045639037, -0.9357086002340218], [1.3760692761371713, 1.9145350254105846, -2.1257143081021397], [0.46117091336276594, 2.8519467789832627, -0.4188289776134376], [2.404864907693512, 3.4780800111829073, -1.0689303403306785], [3.5400592233304673, 1.2295174124845998, -0.43152105569274524], [2.3088468039522576, 1.8960947387583693, 1.2310220414904773], [1.9974224573334862, -0.6906780683055279, 1.0535722235492961], [1.9277183224308914, -0.6529932317206222, -1.1102241252095353], [-0.3501493572534787, 0.5705349971623109, -1.1530217920585815], [-0.4266843221927561, 0.757615307331303, 1.004983428312728], [-4.223062844130574, -2.669723320387231, 2.2813443112863707]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0033', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
