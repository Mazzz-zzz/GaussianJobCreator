import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0067'
logfile = 'conf/5009017845242299296281_0067.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586385, 0.771820394576381, 1.1635336229088504], [-0.3466020415139042, 0.292621294684338, 2.6259377241923803], [-0.9873929842445245, 1.0841115196807707, 3.8194659122849366], [-1.027435612376104, 0.22565869267366792, 5.131350240546323], [-1.1749864386042272, 1.015580149661642, 6.180518858332038], [-2.033642009065594, -0.6428898323589912, 5.085828954561045], [0.5505558709515905, -0.7624012665653352, 5.3666027147286], [1.6529980492054575, 0.07781192399418321, 5.048649815973669], [0.41632863418113125, -1.4669737436897288, 6.58763401860077], [0.3487395885195963, -1.7956834286035943, 4.195218579797673], [-0.2614927667443652, 2.1741733763092714, 4.051469307654177], [-2.230877915023897, 1.4310362312190772, 3.513168225487143], [-0.7784534903451411, -0.9699349115591136, 2.6913872387857145], [0.9731582247379432, 0.30374941440734304, 2.7687950826457888], [-2.0119541879597183, 0.605945531805916, 1.047594832227996], [-0.40994706586378055, 2.0654783699937846, 1.0587059160250767], [1.5770424436171635, 0.0, 0.0], [2.292718146893913, 1.3915527243580577, 0.0], [3.7823355744197107, 1.3186147352454587, -0.4807937764724474], [4.419592497958774, 0.3492507324843974, 0.14976318527135218], [3.8473329977166353, 1.1104432852891981, -1.7795740852228767], [4.368010411051577, 2.468434085335207, -0.20560554420872013], [2.293896732220213, 1.870418904473606, 1.2405689893126328], [1.6494649440008777, 2.235214894314336, -0.8029305726284879], [1.9974224573334805, -0.6906780683055271, 1.0535722235492981], [1.9277183224308942, -0.6529932317206282, -1.1102241252095315], [-0.3501493572534785, -1.283813661620945, 0.08241309473864852], [-0.4266843221927631, 0.49153352563554753, -1.1586058166012272], [0.6931982976311327, -2.6632846875569287, 4.449341856724047]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0067', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
