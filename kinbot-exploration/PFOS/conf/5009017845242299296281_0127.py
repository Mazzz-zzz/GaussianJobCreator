import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0127'
logfile = 'conf/5009017845242299296281_0127.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863847, 0.6217394783082113, -1.2501828803165027], [-0.3976197158559566, -0.07566485901595997, -2.633910119820668], [1.0879414097563251, -0.08367235852192308, -3.1385123060458486], [2.1006320219594947, -0.3015263760374458, -1.960979924321285], [1.6317864231972183, -1.2147664367977424, -1.1287069288319622], [3.280492551617343, -0.6953768716368589, -2.4318418554907217], [2.3884877116331324, 1.276961406080902, -0.9877488342195521], [3.264157839292738, 2.1088893133416184, -1.738434382538826], [1.1272830289191942, 1.692381488956887, -0.49566583848433626], [3.190351489685103, 0.6577071545326659, 0.21810409346137405], [1.2395907752722293, -1.0741885454775733, -4.013081746298284], [1.3717373297746491, 1.0745873687029905, -3.7201686320312395], [-0.7843720924895158, -1.3445493546988834, -2.475221044816154], [-1.1581983763610537, 0.4930222734242895, -3.56149530020636], [-0.25604457595343244, 1.8735740976390407, -1.3916701657561954], [-2.0076024771874503, 0.6427130616946727, -1.0543092166280654], [1.5770424436171646, 0.0, 0.0], [2.292718146893912, 1.3915527243580592, 0.0], [2.3410798567223123, 2.05985269289493, 1.416502376706476], [1.1453916410070282, 2.0331725558915883, 1.9759511228307938], [3.2022179922846403, 1.4422259961598294, 2.1984030628363285], [2.7215555591494893, 3.315205977928799, 1.2745358845394066], [1.6292323391939698, 2.2122558673105788, -0.8090479336198857], [3.5455868300943765, 1.2600392214310814, -0.42809146886198035], [1.997422457333484, -0.6906780683055215, 1.0535722235493004], [1.9277183224308934, -0.6529932317206257, -1.1102241252095306], [-0.3501493572534773, 0.7132786644586332, 1.070608697319931], [-0.4266843221927561, -1.2491488329668505, 0.15362238828850347], [4.142267527937786, 0.73990303960742, 0.06624205314155797]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0127', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
