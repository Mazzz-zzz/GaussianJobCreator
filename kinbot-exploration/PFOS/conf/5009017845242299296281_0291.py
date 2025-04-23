import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0291'
logfile = 'conf/5009017845242299296281_0291.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586379, -1.3935598728846035, 0.08664925740764867], [-2.2709622836291894, -1.3932545648232224, 0.05367636867327392], [-2.9699917885127114, -0.7382819684025499, -1.1888806696804026], [-4.434530032795166, -1.2666412795612016, -1.377425964896209], [-5.118760473256557, -0.4310367516186253, -2.1390661163489013], [-4.423212373540428, -2.4713913694372414, -1.9407634943451875], [-5.334436835556408, -1.426961554161535, 0.2615153052569085], [-4.871160369452701, -2.6048185177058474, 0.9099433253135969], [-5.336387741389749, -0.14210895356265008, 0.8569200082520481], [-6.781729903715371, -1.694958041847532, -0.2988981914185226], [-3.0171454311974566, 0.5789431357488252, -1.0108747386529287], [-2.2839557356380076, -1.0170619338342441, -2.289819390281601], [-2.621205471792921, -2.681987940783288, 0.09005605494849676], [-2.7133779807411256, -0.7960712600074847, 1.1536489463716046], [-0.3710451618282774, -1.9490799840121218, 1.2553873021032798], [-0.24552532002048524, -2.1655580448419816, -0.8973245689258833], [1.577042443617165, 0.0, 0.0], [2.2927181468939164, 1.3915527243580563, 0.0], [3.782335574419711, 1.3186147352454747, -0.48079377647243526], [4.419592497958774, 0.3492507324844033, 0.14976318527135046], [3.8473329977166295, 1.1104432852892154, -1.7795740852228714], [4.368010411051577, 2.468434085335204, -0.20560554420871247], [2.2938967322202166, 1.8704189044736041, 1.2405689893126368], [1.6494649440008744, 2.2352148943143413, -0.8029305726284737], [1.9974224573334853, -0.6906780683055255, 1.0535722235492946], [1.9277183224308951, -0.6529932317206152, -1.1102241252095357], [-0.3501493572534774, 0.5705349971623099, -1.153021792058573], [-0.4266843221927577, 0.7576153073312948, 1.004983428312732], [-6.958333004020121, -2.645553984935081, -0.3326781358271219]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0291', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
