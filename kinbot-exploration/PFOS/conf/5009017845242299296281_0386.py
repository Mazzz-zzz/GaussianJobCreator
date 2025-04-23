import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0386'
logfile = 'conf/5009017845242299296281_0386.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863827, 0.6217394783082154, -1.250182880316501], [-2.2709622836291934, 0.6501421835576524, -1.2334320314121736], [-3.0203184489305848, -0.727198632127044, -1.1788705129599804], [-4.480058160568046, -0.5732592173009051, -0.626559196912983], [-4.453761713276069, -0.47762555477201196, 0.6912665284878111], [-5.05670016854204, 0.5106696876177764, -1.1379124645025485], [-5.555811090753356, -2.042513255532673, -1.0801135326351357], [-6.700610364745075, -2.032842561010063, -0.23644765816498706], [-5.612894691041407, -2.0893298783299, -2.4942940154788342], [-4.587105203609392, -3.189372367441306, -0.6038093110983515], [-3.086116569243987, -1.2298413735800906, -2.4085096609760925], [-2.36392887439246, -1.5654026010903253, -0.3870393224017533], [-2.5961406643712652, 1.3384499968011863, -0.135561883504018], [-2.6897770261666536, 1.3195442638618289, -2.300575080214779], [-0.3710451618282835, -0.1126573032037982, -2.3156464312138967], [-0.24552532002049293, 1.8598848945507231, -1.4267659957399772], [1.5770424436171668, 0.0, 0.0], [2.292718146893916, 1.39155272435806, 0.0], [2.3410798567223208, 2.059852692894929, 1.4165023767064773], [1.1453916410070362, 2.033172555891591, 1.9759511228307938], [3.202217992284648, 1.4422259961598325, 2.198403062836328], [2.7215555591494978, 3.315205977928801, 1.2745358845394033], [1.6292323391939723, 2.2122558673105805, -0.8090479336198845], [3.5455868300943787, 1.2600392214310843, -0.42809146886198585], [1.997422457333486, -0.6906780683055243, 1.053572223549291], [1.9277183224308951, -0.6529932317206222, -1.1102241252095346], [-0.3501493572534733, 0.7132786644586351, 1.0706086973199336], [-0.4266843221927547, -1.249148832966849, 0.1536223882884999], [-3.6713070194980872, -2.8774814265258564, -0.6041464302345645]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0386', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
