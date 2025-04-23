import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0365'
logfile = 'conf/5009017845242299296281_0365.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863847, 0.6217394783082113, -1.2501828803165027], [-0.3976197158559566, -0.07566485901595997, -2.633910119820668], [1.0879414097563251, -0.08367235852192308, -3.1385123060458486], [2.1006320219594947, -0.3015263760374458, -1.960979924321285], [1.6317864231972183, -1.2147664367977424, -1.1287069288319622], [3.280492551617343, -0.6953768716368589, -2.4318418554907217], [2.3884877116331324, 1.276961406080902, -0.9877488342195521], [3.264157839292738, 2.1088893133416193, -1.738434382538825], [1.1272830289191944, 1.6923814889568873, -0.4956658384843364], [3.190351489685103, 0.6577071545326659, 0.21810409346137405], [1.2395907752722293, -1.0741885454775733, -4.013081746298284], [1.3717373297746491, 1.0745873687029905, -3.7201686320312395], [-0.7843720924895158, -1.3445493546988834, -2.475221044816154], [-1.1581983763610537, 0.4930222734242895, -3.56149530020636], [-0.25604457595343244, 1.8735740976390407, -1.3916701657561954], [-2.0076024771874503, 0.6427130616946727, -1.0543092166280654], [1.5770424436171646, 0.0, 0.0], [2.292718146893912, 1.3915527243580592, 0.0], [3.782335574419712, 1.3186147352454631, -0.48079377647244265], [4.419592497958772, 0.34925073248439675, 0.14976318527135185], [3.8473329977166335, 1.1104432852892019, -1.7795740852228785], [4.368010411051574, 2.468434085335211, -0.20560554420871807], [2.293896732220212, 1.87041890447361, 1.240568989312632], [1.649464944000876, 2.2352148943143373, -0.8029305726284884], [1.9974224573334838, -0.6906780683055213, 1.0535722235493004], [1.9277183224308938, -0.652993231720626, -1.1102241252095302], [-0.3501493572534773, 0.7132786644586332, 1.070608697319931], [-0.4266843221927561, -1.2491488329668505, 0.15362238828850347], [2.9682483780573006, 1.1165443981350174, 1.040357690829392]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0365', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
