import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0055'
logfile = 'conf/5009017845242299296281_0055.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863802, 0.6217394783082122, -1.2501828803165043], [-0.3976197158559543, -0.07566485901596291, -2.6339101198206682], [-0.7491833517666566, -1.5988198371506752, -2.768004365822093], [0.37513150925004174, -2.51023981580582, -2.164005297273753], [1.3770046036719104, -2.6139538171368173, -3.019572741534261], [0.8207913513067474, -1.9990413104241207, -1.019905898695764], [-0.25282578633505265, -4.239487159013335, -1.7938207594467366], [-1.0611560354132707, -4.662186520076783, -2.8848562613400928], [0.8456628214350043, -4.976343107004426, -1.288130278858511], [-1.1898482641172174, -3.8858833030415476, -0.5782979011746819], [-1.8820373528306358, -1.8378201814607742, -2.113395933093805], [-0.8963077255402084, -1.9155699118782956, -4.048025307170628], [-1.1305489784219294, 0.5952337900285942, -3.526959494450472], [0.8878652548597795, 0.09037708117405727, -2.9207161247297275], [-0.25604457595342645, 1.8735740976390398, -1.3916701657562012], [-2.0076024771874446, 0.6427130616946735, -1.0543092166280668], [1.5770424436171646, 0.0, 0.0], [2.2927181468939137, 1.3915527243580552, 0.0], [2.3410798567223203, 2.0598526928949266, 1.4165023767064762], [1.1453916410070306, 2.033172555891592, 1.9759511228307929], [3.202217992284642, 1.4422259961598263, 2.198403062836328], [2.721555559149496, 3.315205977928796, 1.2745358845394041], [1.6292323391939756, 2.2122558673105797, -0.8090479336198867], [3.5455868300943796, 1.2600392214310747, -0.4280914688619768], [1.9974224573334822, -0.6906780683055253, 1.0535722235493021], [1.927718322430895, -0.6529932317206275, -1.110224125209532], [-0.3501493572534767, 0.7132786644586364, 1.0706086973199318], [-0.4266843221927588, -1.2491488329668496, 0.15362238828850447], [-1.1717848652322163, -4.596617931856037, 0.07782128595571074]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0055', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
