import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0488'
logfile = 'conf/5009017845242299296281_0488.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863837, 0.6217394783082166, -1.2501828803164992], [-0.3976197158559577, -0.07566485901595234, -2.633910119820667], [-1.1233882121466363, 0.4940238789133735, -3.902883790150536], [-0.3678111583333681, 0.11496132015470231, -5.223868290930236], [0.05889782664840387, -1.1338827048026767, -5.154514224423079], [-1.1736420480321825, 0.2500793113542995, -6.273287056826106], [1.1196427974117333, 1.220622106908921, -5.518246084398577], [1.9616071617376694, 0.5883203812387445, -6.474133915019792], [0.6301620320434117, 2.543463970565568, -5.644161753557944], [1.7841101462045712, 1.0914472580939258, -4.096156203051111], [-2.351535285418468, -0.013330338437761258, -3.959787366978748], [-1.18853011037124, 1.8171105205812643, -3.8278633047836688], [0.9186702629217393, 0.05926285465343883, -2.818308304723229], [-0.6863393184328443, -1.3666682794365088, -2.5232286857846637], [-0.2560445759534295, 1.8735740976390445, -1.3916701657561927], [-2.0076024771874468, 0.6427130616946789, -1.054309216628057], [1.5770424436171657, 0.0, 0.0], [2.2927181468939164, 1.3915527243580565, 0.0], [3.782335574419713, 1.3186147352454654, -0.4807937764724456], [4.419592497958778, 0.3492507324843952, 0.14976318527134957], [3.847332997716634, 1.1104432852892008, -1.7795740852228787], [4.3680104110515785, 2.4684340853352094, -0.2056055442087198], [2.2938967322202175, 1.8704189044736066, 1.2405689893126353], [1.6494649440008797, 2.2352148943143417, -0.8029305726284853], [1.9974224573334838, -0.690678068305525, 1.0535722235492988], [1.927718322430894, -0.6529932317206194, -1.1102241252095326], [-0.3501493572534726, 0.7132786644586362, 1.0706086973199354], [-0.4266843221927547, -1.249148832966851, 0.15362238828850353], [2.1940624879620567, 1.9292923024268884, -3.839412856749302]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0488', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
