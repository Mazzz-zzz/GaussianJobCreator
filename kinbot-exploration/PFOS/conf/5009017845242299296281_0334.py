import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0334'
logfile = 'conf/5009017845242299296281_0334.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586382, 0.771820394576386, 1.1635336229088489], [-2.270962283629193, 0.7431123812655659, 1.1797556627388994], [-3.020318448930583, 1.3845311280592616, -0.040337232539323836], [-2.282344703997105, 1.094621548070399, -1.3934567195845522], [-3.1133589271330373, 1.2733031875488667, -2.405380716600173], [-1.2352156515326023, 1.9028124788375202, -1.53219462883643], [-1.6404353738143138, -0.6669399139804485, -1.4741843998222144], [-1.3547325602404028, -0.9743817283410166, -2.8329648975144037], [-0.7121068992865005, -0.8152132836813629, -0.4151385920941372], [-2.964071167648961, -1.4107992895555208, -1.0556104583760326], [-4.245418844080689, 0.8714978469683302, -0.1104615992182993], [-3.0940580938002604, 2.6998070798223184, 0.11790763144453659], [-2.5961406643712635, -0.5518249635012413, 1.2269126406770399], [-2.6897770261666487, 1.332584330848512, 2.2930463940297705], [-0.37104516182828245, 2.061737287215913, 1.0602591291106098], [-0.2455253200204891, 0.3056731502912654, 2.324090564665854], [1.5770424436171668, 0.0, 0.0], [2.292718146893917, 1.3915527243580572, 0.0], [2.3410798567223208, 2.0598526928949297, 1.4165023767064746], [1.1453916410070364, 2.033172555891593, 1.9759511228307909], [3.202217992284643, 1.4422259961598294, 2.198403062836328], [2.721555559149501, 3.3152059779287972, 1.2745358845394041], [1.6292323391939774, 2.212255867310575, -0.8090479336198866], [3.5455868300943827, 1.2600392214310774, -0.42809146886198085], [1.9974224573334838, -0.6906780683055275, 1.053572223549301], [1.927718322430896, -0.6529932317206257, -1.110224125209532], [-0.3501493572534751, -1.2838136616209437, 0.08241309473864959], [-0.4266843221927562, 0.4915335256355431, -1.1586058166012267], [-3.451078465921009, -1.6975845533152814, -1.8408113820993843]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0334', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
