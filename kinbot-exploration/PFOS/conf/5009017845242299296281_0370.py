import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0370'
logfile = 'conf/5009017845242299296281_0370.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863863, 0.6217394783082137, -1.2501828803165034], [-0.3976197158559574, -0.07566485901595767, -2.633910119820669], [1.0879414097563262, -0.0836723585219217, -3.1385123060458495], [1.473042089945535, 1.268447375297704, -3.833533780919367], [2.7879690257787106, 1.3947914162119797, -3.8723446018695187], [0.9875700214111528, 1.3076225145166924, -5.071145278802743], [0.7750378848301898, 2.7483816085676347, -2.9147040219585287], [-0.6091697343464533, 2.8516209656536944, -3.2241259183585473], [1.2916422504512648, 2.69546051552042, -1.5972534312417703], [1.5508075901178944, 3.881525072116208, -3.685878376552856], [1.8947536492841406, -0.2666142022536076, -2.0970816341016425], [1.2598103705638566, -1.064159071597839, -4.015744696739172], [-0.7843720924895178, -1.3445493546988798, -2.4752210448161556], [-1.1581983763610542, 0.49302227342429, -3.561495300206358], [-0.2560445759534314, 1.873574097639041, -1.391670165756194], [-2.0076024771874486, 0.6427130616946721, -1.0543092166280648], [1.5770424436171633, 0.0, 0.0], [2.2927181468939124, 1.3915527243580568, 0.0], [1.6005215470082481, 2.4407219045638993, -0.9357086002340306], [1.376069276137161, 1.9145350254105782, -2.1257143081021446], [0.4611709133627515, 2.85194677898325, -0.4188289776134486], [2.404864907693499, 3.478080011182909, -1.068930340330683], [3.540059223330463, 1.2295174124846113, -0.43152105569274574], [2.3088468039522474, 1.896094738758375, 1.2310220414904693], [1.997422457333481, -0.6906780683055215, 1.0535722235493], [1.9277183224308945, -0.6529932317206217, -1.110224125209531], [-0.3501493572534754, 0.7132786644586363, 1.0706086973199316], [-0.4266843221927558, -1.2491488329668532, 0.15362238828850042], [2.365831876474917, 3.529203820234931, -4.070023576727266]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0370', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
