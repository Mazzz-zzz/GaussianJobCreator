import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0035'
logfile = 'conf/5009017845242299296281_0035.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863756, -1.3935598728846026, 0.08664925740765234], [-2.2709622836291863, -1.3932545648232288, 0.05367636867327788], [-3.0203184489305785, -0.6573324959322318, 1.219207745499312], [-3.131423238500171, -1.5575296086444992, 2.49861621936573], [-4.113036047167041, -2.4298037535513113, 2.3500360959803137], [-1.9898347010197868, -2.2068663000332647, 2.708467387727265], [-3.485344517394217, -0.5425193845977196, 4.036848095430425], [-2.275694716609433, 0.08480777878383478, 4.443835806983219], [-4.69823536319325, 0.14983481895219408, 3.8025443936848498], [-3.796615807083963, -1.7284975645804368, 5.025336117931851], [-2.3454389866601097, 0.4439741098051188, 1.5364663530591403], [-4.247573342377318, -0.3344495623799493, 0.8318812175920427], [-2.5961406643712612, -0.7866250332999511, -1.0913507571730114], [-2.689777026166637, -2.652128594710351, 0.007528686185015135], [-0.3710451618282712, -1.949079984012114, 1.2553873021032869], [-0.24552532002048144, -2.1655580448419838, -0.8973245689258763], [1.5770424436171677, 0.0, 0.0], [2.2927181468939146, 1.3915527243580597, 0.0], [2.3410798567223248, 2.059852692894932, 1.4165023767064726], [1.1453916410070417, 2.0331725558915963, 1.9759511228307929], [3.2022179922846554, 1.4422259961598336, 2.1984030628363236], [2.7215555591494947, 3.3152059779288052, 1.2745358845394], [1.629232339193969, 2.2122558673105788, -0.8090479336198849], [3.5455868300943827, 1.2600392214310836, -0.42809146886198846], [1.9974224573334882, -0.6906780683055223, 1.0535722235492964], [1.9277183224308962, -0.652993231720618, -1.1102241252095357], [-0.35014935725348056, 0.5705349971623069, -1.1530217920585812], [-0.4266843221927565, 0.7576153073312993, 1.0049834283127286], [-4.493119589630578, -1.4740712987633844, 5.646715344215326]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0035', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
