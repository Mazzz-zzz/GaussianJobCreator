import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0378'
logfile = 'conf/5009017845242299296281_0378.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863828, -1.3935598728845986, 0.08664925740765334], [-0.3466020415139064, -2.4204394252486647, -1.0595513872112596], [-0.7363023803695227, -2.0273092602547935, -2.527598274606108], [0.11232051632788868, -2.8153848177000156, -3.5851242808869475], [0.2504583136006516, -4.072929646905836, -3.2030236467873316], [-0.4794659588750192, -2.776614081266731, -4.775580344047762], [1.8265739031510893, -2.0802541831649792, -3.7914773934108306], [2.3338654001233095, -1.7654983828297763, -2.500735307994452], [2.5077972431396933, -2.88463741650845, -4.737130155590843], [1.4148838370527916, -0.7418952538065334, -4.512483146155908], [-2.0209682759934697, -2.311513884043145, -2.7221562889464304], [-0.5312999584009223, -0.7299509739637069, -2.715388810757052], [0.9809271994195052, -2.562327183726283, -1.010956477982444], [-0.9145558447621731, -3.5839764719224805, -0.7664486088494572], [-2.0119541879597196, -1.2102165034856958, 0.0009668077396009832], [-0.4099470658637758, -1.9496054034114851, 1.2594037813693548], [1.5770424436171633, 0.0, 0.0], [2.2927181468939146, 1.3915527243580539, 0.0], [3.7823355744197187, 1.318614735245457, -0.48079377647244337], [4.419592497958775, 0.34925073248439176, 0.1497631852713457], [3.8473329977166344, 1.1104432852891992, -1.7795740852228805], [4.368010411051581, 2.4684340853351987, -0.20560554420871946], [2.2938967322202175, 1.8704189044736004, 1.2405689893126344], [1.6494649440008802, 2.2352148943143373, -0.8029305726284827], [1.997422457333483, -0.6906780683055285, 1.0535722235492984], [1.927718322430894, -0.6529932317206261, -1.110224125209534], [-0.3501493572534776, 0.5705349971623103, -1.15302179205858], [-0.42668432219275776, 0.7576153073313008, 1.0049834283127324], [1.3340244013882485, -0.024402921307010057, -3.868566310559571]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0378', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
