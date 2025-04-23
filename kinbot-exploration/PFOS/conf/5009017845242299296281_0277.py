import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0277'
logfile = 'conf/5009017845242299296281_0277.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863837, 0.6217394783082144, -1.2501828803165003], [-2.2709622836291934, 0.6501421835576524, -1.2334320314121736], [-2.9699917885127123, 1.3987418462127572, -0.04493060495238103], [-2.195640795695414, 2.7019843977942153, 0.356970565327479], [-1.781466718272871, 3.331102206613142, -0.7289265565468985], [-2.982357122914438, 3.5087863141483293, 1.0634039512934532], [-0.6934923395346748, 2.3221193673769, 1.4157037259656842], [-0.03155894258452733, 1.192880523077734, 0.8597407610679074], [-0.0593969656800988, 3.5585968531411893, 1.688473327831917], [-1.445445794483658, 1.8729277458066866, 2.7246481550803843], [-4.201534988011382, 1.7379422988872675, -0.4153843175535463], [-3.0244357979279344, 0.6050007412875026, 1.0168795050389612], [-2.621205471792929, 1.2630031390416294, -2.367697716836086], [-2.71337798074113, -0.6010536646032236, -1.2662424075749565], [-0.37104516182828257, -0.11265730320379891, -2.3156464312139007], [-0.24552532002049013, 1.85988489455072, -1.4267659957399774], [1.5770424436171642, 0.0, 0.0], [2.292718146893917, 1.3915527243580588, 0.0], [2.3410798567223248, 2.059852692894926, 1.4165023767064777], [1.145391641007037, 2.0331725558915945, 1.9759511228307958], [3.202217992284647, 1.4422259961598287, 2.1984030628363262], [2.721555559149503, 3.3152059779288, 1.2745358845394033], [1.6292323391939747, 2.212255867310578, -0.8090479336198837], [3.545586830094381, 1.260039221431077, -0.4280914688619833], [1.9974224573334827, -0.6906780683055288, 1.0535722235492997], [1.9277183224308956, -0.6529932317206222, -1.1102241252095348], [-0.35014935725347535, 0.7132786644586341, 1.0706086973199345], [-0.42668432219275787, -1.2491488329668514, 0.1536223882885013], [-1.5345648316065033, 0.9098563549747675, 2.747312666532669]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0277', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
