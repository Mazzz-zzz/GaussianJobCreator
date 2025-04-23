import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0394'
logfile = 'conf/5009017845242299296281_0394.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863827, 0.6217394783082104, -1.2501828803165036], [-0.3976197158559555, -0.0756648590159635, -2.6339101198206682], [-0.7491833517666558, -1.5988198371506765, -2.768004365822093], [-2.0819874597153203, -1.9555730430184177, -2.022386823820662], [-1.8542445281745772, -2.0700904867547636, -0.7256521341904679], [-2.9956443258697933, -1.011627170194836, -2.229985862019596], [-2.8081615356386624, -3.5761500072006207, -2.629070304503851], [-1.7481451191270514, -4.513827780983398, -2.7689152081994903], [-3.9802865149534536, -3.8112279263822417, -1.869979943329676], [-3.2557013693023076, -3.102288475544046, -4.062843639417499], [-0.895313736421605, -1.8945036554924028, -4.056495166238397], [0.2259041467951681, -2.333069003291668, -2.2478963639044998], [-1.1305489784219298, 0.5952337900285927, -3.5269594944504723], [0.8878652548597811, 0.09037708117405664, -2.9207161247297284], [-0.2560445759534277, 1.8735740976390403, -1.391670165756204], [-2.0076024771874463, 0.64271306169467, -1.0543092166280654], [1.5770424436171662, 0.0, 0.0], [2.292718146893918, 1.3915527243580537, 0.0], [2.3410798567223323, 2.0598526928949283, 1.4165023767064704], [1.1453916410070433, 2.0331725558915985, 1.9759511228307884], [3.2022179922846528, 1.4422259961598232, 2.1984030628363236], [2.7215555591495173, 3.3152059779287937, 1.274535884539396], [1.6292323391939871, 2.2122558673105774, -0.8090479336198891], [3.545586830094387, 1.2600392214310667, -0.4280914688619808], [1.9974224573334853, -0.6906780683055251, 1.0535722235493006], [1.9277183224308987, -0.6529932317206296, -1.1102241252095302], [-0.35014935725347457, 0.7132786644586374, 1.0706086973199311], [-0.4266843221927564, -1.2491488329668508, 0.1536223882885046], [-2.5761887559051284, -3.3201626257970376, -4.716109167041354]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0394', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
