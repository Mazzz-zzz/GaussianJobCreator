import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0298'
logfile = 'conf/5009017845242299296281_0298.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863816, 0.7718203945763855, 1.1635336229088482], [-0.3466020415139031, 0.2926212946843427, 2.6259377241923807], [1.1624535490467065, 0.30841610108309514, 3.0549844500983023], [1.9243606525366, -0.9650386511798187, 2.547642633772617], [3.226844590515571, -0.741451815104614, 2.5544588730337825], [1.6500667572714027, -2.008979336510665, 3.324720972863477], [1.4242023469536254, -1.4263133994933446, 0.7987199347834181], [1.3520611618988883, -0.23362569436472597, 0.027537274667205965], [2.192125523486113, -2.559765572589417, 0.4368678408333378], [-0.04001359381904091, -1.9199612818331857, 1.103557722561794], [1.7462981079693127, 1.383924747436088, 2.5340254928450245], [1.2535978601361863, 0.34611839988490517, 4.378125052116878], [-1.0117971004258957, 1.1239736000383738, 3.4329051118380085], [-0.8215282967123796, -0.9363359387499647, 2.7884443111439534], [-2.011954187959717, 0.6059455318059157, 1.047594832227993], [-0.40994706586377827, 2.0654783699937846, 1.0587059160250734], [1.5770424436171666, 0.0, 0.0], [2.2927181468939146, 1.3915527243580579, 0.0], [1.6005215470082463, 2.440721904563902, -0.9357086002340322], [1.3760692761371618, 1.9145350254105735, -2.125714308102147], [0.4611709133627484, 2.851946778983252, -0.4188289776134525], [2.4048649076934905, 3.4780800111829135, -1.0689303403306907], [3.5400592233304646, 1.2295174124846127, -0.4315210556927457], [2.3088468039522483, 1.8960947387583786, 1.2310220414904671], [1.9974224573334847, -0.6906780683055236, 1.0535722235492995], [1.927718322430899, -0.6529932317206233, -1.1102241252095313], [-0.3501493572534738, -1.2838136616209432, 0.08241309473864956], [-0.4266843221927561, 0.4915335256355466, -1.1586058166012279], [-0.1134427660932617, -2.205310393380818, 2.0250491375484434]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0298', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
