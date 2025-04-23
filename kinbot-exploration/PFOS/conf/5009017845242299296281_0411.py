import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0411'
logfile = 'conf/5009017845242299296281_0411.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.69372834458638, -1.3935598728846024, 0.08664925740765357], [-2.27096228362919, -1.3932545648232229, 0.05367636867327765], [-2.9699917885127127, -0.7382819684025512, -1.1888806696803988], [-4.434530032795169, -1.2666412795612028, -1.3774259648962028], [-4.4121494903032055, -2.4560673876137655, -1.952983597004151], [-5.0498579838776125, -1.3550224808038638, -0.20164523855451022], [-5.456246996148601, -0.11973325061456805, -2.455553856682976], [-6.600925948474646, -0.8340666475877683, -2.904852089353076], [-5.511129515012079, 1.1292950530658474, -1.7905177844459568], [-4.446822673886868, 0.023241914471068055, -3.6560384545034434], [-3.0171454311974584, 0.5789431357488226, -1.0108747386529275], [-2.283955735638008, -1.0170619338342468, -2.2898193902815973], [-2.621205471792919, -2.6819879407832863, 0.09005605494850148], [-2.713377980741126, -0.7960712600074806, 1.1536489463716078], [-0.37104516182827785, -1.949079984012116, 1.255387302103288], [-0.2455253200204861, -2.165558044841982, -0.8973245689258778], [1.577042443617164, 0.0, 0.0], [2.2927181468939137, 1.3915527243580585, 0.0], [2.3410798567223243, 2.05985269289494, 1.4165023767064648], [1.145391641007037, 2.033172555891596, 1.975951122830791], [3.202217992284651, 1.442225996159835, 2.1984030628363165], [2.72155555914949, 3.3152059779287946, 1.2745358845394037], [1.6292323391939698, 2.2122558673105766, -0.8090479336198838], [3.5455868300943743, 1.2600392214310803, -0.4280914688619907], [1.997422457333486, -0.6906780683055195, 1.0535722235492975], [1.927718322430894, -0.652993231720618, -1.110224125209533], [-0.35014935725347746, 0.570534997162309, -1.1530217920585724], [-0.42668432219275865, 0.7576153073312981, 1.0049834283127308], [-3.5429136097678424, -0.14920133061214555, -3.3574180973195467]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0411', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
