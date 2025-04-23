import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0409'
logfile = 'conf/5009017845242299296281_0409.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863799, 0.7718203945763888, 1.1635336229088442], [-0.3466020415139041, 0.2926212946843473, 2.625937724192379], [1.162453549046706, 0.3084161010830979, 3.0549844500983028], [1.9172819807099775, 1.5569601130304864, 2.479695785178602], [2.2272170921886345, 1.3494186210778603, 1.2118753330372962], [1.1566135117755332, 2.643560915797122, 2.577584581007528], [3.50626027001411, 1.903654662880097, 3.4159725692225544], [4.313987042727912, 2.761027758671974, 2.619153993219461], [3.1421067329679917, 2.1770371113498426, 4.756871898021516], [4.116322039547438, 0.45223870139642713, 3.374032446269606], [1.2315520144708516, 0.3493659548263257, 4.382591343745828], [1.768541893208141, -0.7841917084166116, 2.6085598817993625], [-1.0117971004258957, 1.1239736000383767, 3.432905111838007], [-0.8215282967123818, -0.9363359387499599, 2.788444311143955], [-2.011954187959715, 0.6059455318059207, 1.0475948322279887], [-0.40994706586377555, 2.065478369993788, 1.0587059160250727], [1.5770424436171662, 0.0, 0.0], [2.2927181468939177, 1.3915527243580543, 0.0], [3.782335574419716, 1.3186147352454596, -0.48079377647244925], [4.419592497958776, 0.34925073248438765, 0.14976318527135069], [3.8473329977166384, 1.1104432852891963, -1.7795740852228783], [4.36801041105158, 2.468434085335204, -0.20560554420871946], [2.293896732220218, 1.870418904473607, 1.2405689893126308], [1.649464944000888, 2.235214894314334, -0.8029305726284889], [1.9974224573334838, -0.6906780683055276, 1.0535722235493], [1.9277183224308982, -0.652993231720632, -1.1102241252095293], [-0.3501493572534737, -1.2838136616209452, 0.08241309473865081], [-0.426684322192755, 0.4915335256355452, -1.158605816601228], [3.4139916343943835, -0.2045601047302396, 3.267685293208864]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0409', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
