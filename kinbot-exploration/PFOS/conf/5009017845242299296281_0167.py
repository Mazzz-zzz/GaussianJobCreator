import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0167'
logfile = 'conf/5009017845242299296281_0167.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863823, 0.7718203945763816, 1.1635336229088493], [-0.34660204151390084, 0.2926212946843405, 2.6259377241923816], [-0.9873929842445218, 1.084111519680774, 3.819465912284936], [-0.1847769512407465, 2.390826822806651, 4.147928663816883], [0.8897635883676669, 2.0915643040472673, 4.856664522004434], [0.18128897051813278, 3.0037758903973177, 3.0257926551175993], [-1.2119219398679448, 3.609697700936593, 5.13837820109769], [-2.1177150707628143, 4.252726397317443, 4.2504120376295], [-1.6023819083198627, 2.9473596998896188, 6.327582202543041], [-0.06536911483788746, 4.623001131115862, 5.511440404505566], [-2.2305068985543333, 1.4246952716906078, 3.4913469454112684], [-1.003631458776227, 0.3235218506680348, 4.906512125113522], [-0.7784534903451336, -0.9699349115591119, 2.691387238785717], [0.9731582247379497, 0.3037494144073529, 2.768795082645788], [-2.0119541879597156, 0.6059455318059139, 1.0475948322279944], [-0.40994706586378016, 2.0654783699937824, 1.0587059160250756], [1.5770424436171675, 0.0, 0.0], [2.292718146893914, 1.3915527243580594, 0.0], [2.3410798567223177, 2.0598526928949266, 1.4165023767064764], [1.1453916410070342, 2.0331725558915905, 1.975951122830794], [3.202217992284642, 1.442225996159831, 2.19840306283633], [2.7215555591494978, 3.315205977928797, 1.2745358845394066], [1.629232339193969, 2.2122558673105783, -0.8090479336198833], [3.54558683009438, 1.2600392214310812, -0.4280914688619817], [1.9974224573334864, -0.690678068305524, 1.0535722235493], [1.9277183224308996, -0.6529932317206261, -1.1102241252095295], [-0.3501493572534705, -1.2838136616209446, 0.08241309473865077], [-0.4266843221927554, 0.49153352563554276, -1.1586058166012276], [-0.048724513182227744, 5.358017802990644, 4.8826035130958525]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0167', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
