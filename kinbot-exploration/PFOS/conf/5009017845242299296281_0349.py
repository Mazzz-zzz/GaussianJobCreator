import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0349'
logfile = 'conf/5009017845242299296281_0349.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863798, 0.7718203945763887, 1.163533622908844], [-0.34660204151390284, 0.29262129468434483, 2.6259377241923785], [1.1624535490467067, 0.30841610108309786, 3.0549844500983023], [1.9172819807099764, 1.556960113030487, 2.4796957851786035], [2.2272170921886336, 1.349418621077861, 1.2118753330372967], [1.1566135117755314, 2.6435609157971225, 2.5775845810075286], [3.50626027001411, 1.903654662880097, 3.4159725692225544], [4.165236172296549, 0.664924799205824, 3.6477700947972567], [4.114144233530834, 3.0215131987428414, 2.7945016106524743], [2.8697492907789988, 2.398056333467619, 4.769103527867149], [1.2315520144708525, 0.3493659548263219, 4.382591343745828], [1.768541893208142, -0.7841917084166112, 2.608559881799361], [-1.0117971004258957, 1.1239736000383767, 3.432905111838007], [-0.8215282967123823, -0.9363359387499601, 2.788444311143956], [-2.0119541879597147, 0.6059455318059188, 1.0475948322279898], [-0.409947065863776, 2.0654783699937873, 1.0587059160250736], [1.5770424436171662, 0.0, 0.0], [2.2927181468939177, 1.3915527243580557, 0.0], [1.6005215470082506, 2.4407219045639006, -0.9357086002340272], [1.3760692761371667, 1.914535025410577, -2.1257143081021432], [0.4611709133627562, 2.851946778983253, -0.41882897761344906], [2.4048649076935025, 3.4780800111829104, -1.0689303403306838], [3.5400592233304664, 1.2295174124846064, -0.4315210556927444], [2.3088468039522505, 1.8960947387583755, 1.2310220414904665], [1.9974224573334831, -0.6906780683055264, 1.0535722235492995], [1.9277183224308982, -0.6529932317206302, -1.1102241252095304], [-0.35014935725347357, -1.283813661620946, 0.08241309473865085], [-0.42668432219275537, 0.49153352563554575, -1.158605816601228], [2.8170694699288066, 1.668330363895973, 5.402102131739162]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0349', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
