import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0422'
logfile = 'conf/5009017845242299296281_0422.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863829, 0.7718203945763835, 1.1635336229088495], [-2.2709622836291934, 0.7431123812655637, 1.179755662738898], [-3.0203184489305848, 1.384531128059259, -0.040337232539325606], [-4.480058160568045, 0.8292457901518746, -0.18317744667967567], [-5.037241623014006, 0.7188730088199354, 1.0100812768982363], [-5.212262690872638, 1.6411358795975055, -0.9405485901848986], [-4.5062276883420225, -0.8628257980321968, -0.9942706764574425], [-4.31336332463281, -0.6888607354849567, -2.392473273377493], [-3.71978584464473, -1.7187344516148035, -0.18539097288438816], [-6.018534089092507, -1.2168577934012865, -0.7332470638989184], [-3.086116569243988, 2.7007512384555885, 0.13918095834253655], [-2.363928874392461, 1.1178871860085953, -1.162158758493585], [-2.5961406643712643, -0.551824963501243, 1.2269126406770379], [-2.6897770261666523, 1.3325843308485132, 2.2930463940297683], [-0.3710451618282815, 2.061737287215912, 1.0602591291106067], [-0.24552532002049293, 0.3056731502912636, 2.324090564665856], [1.5770424436171662, 0.0, 0.0], [2.2927181468939173, 1.3915527243580543, 0.0], [2.3410798567223248, 2.059852692894929, 1.4165023767064715], [1.1453916410070393, 2.033172555891594, 1.9759511228307907], [3.202217992284644, 1.442225996159828, 2.1984030628363254], [2.721555559149507, 3.3152059779287977, 1.2745358845394004], [1.6292323391939783, 2.212255867310575, -0.8090479336198921], [3.5455868300943827, 1.2600392214310707, -0.42809146886198357], [1.9974224573334842, -0.6906780683055243, 1.0535722235493024], [1.927718322430898, -0.6529932317206244, -1.1102241252095286], [-0.3501493572534737, -1.2838136616209428, 0.08241309473865066], [-0.42668432219275576, 0.4915335256355448, -1.1586058166012259], [-6.55502152370259, -0.9896134390280463, -1.50558353371436]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0422', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
