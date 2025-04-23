import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0255'
logfile = 'conf/5009017845242299296281_0255.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863814, 0.6217394783082147, -1.250182880316504], [-0.3976197158559531, -0.07566485901595822, -2.6339101198206682], [-0.7491833517666566, -1.59881983715067, -2.768004365822096], [-0.9314389667275241, -2.023554908144312, -4.266659819649047], [-0.00874900226188584, -1.442588222006411, -5.013401329242118], [-0.833173506702976, -3.3442788342087293, -4.3888574098059925], [-2.614049859498324, -1.5265903701332353, -4.9326254063000965], [-2.5614370490999603, -1.5723755682523104, -6.353036127051377], [-3.5817771310013504, -2.204825148613837, -4.152331689231532], [-2.6037140182008875, -0.019548424781951387, -4.47517554550057], [0.24079307197239197, -2.316267212136121, -2.2443561088085517], [-1.8795656402123522, -1.858056375972156, -2.123440869805504], [-1.1305489784219271, 0.5952337900286016, -3.5269594944504723], [0.8878652548597826, 0.09037708117406053, -2.920716124729728], [-0.2560445759534255, 1.873574097639044, -1.3916701657561978], [-2.007602477187444, 0.6427130616946777, -1.0543092166280663], [1.5770424436171657, 0.0, 0.0], [2.2927181468939164, 1.3915527243580548, 0.0], [3.782335574419716, 1.3186147352454571, -0.4807937764724435], [4.419592497958775, 0.34925073248438865, 0.14976318527135368], [3.847332997716639, 1.1104432852891923, -1.7795740852228725], [4.3680104110515785, 2.4684340853352023, -0.2056055442087183], [2.293896732220219, 1.8704189044736077, 1.2405689893126317], [1.6494649440008846, 2.235214894314337, -0.8029305726284901], [1.9974224573334847, -0.6906780683055251, 1.0535722235492984], [1.9277183224308954, -0.6529932317206278, -1.1102241252095295], [-0.35014935725347257, 0.7132786644586405, 1.0706086973199296], [-0.42668432219275804, -1.2491488329668474, 0.15362238828850197], [-2.009639682285006, 0.09860050126156159, -3.720802361077325]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0255', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
