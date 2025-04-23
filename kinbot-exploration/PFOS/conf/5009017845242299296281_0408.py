import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0408'
logfile = 'conf/5009017845242299296281_0408.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863814, 0.6217394783082119, -1.2501828803165036], [-2.270962283629193, 0.6501421835576519, -1.2334320314121778], [-2.969991788512709, 1.398741846212765, -0.04493060495238687], [-4.434530032795165, 1.8262065172130126, -0.4082305431339223], [-5.030583914311105, 0.8584013505576266, -1.0824522224070845], [-5.127712885957733, 2.0886527649130584, 0.6960880430098529], [-4.4665010348623895, 3.375387929270162, -1.4668514094871372], [-3.447737270510894, 3.261618563650534, -2.452553155301199], [-5.825569276589183, 3.6463126376074753, -1.7582104463986359], [-3.9983000348796893, 4.407779536519637, -0.373463087478414], [-3.017145431197455, 0.5859716358429901, 1.006816832231575], [-2.2839557356380036, 2.4915727289791874, 0.2641082232182081], [-2.6212054717929245, 1.2630031390416254, -2.3676977168360884], [-2.7133779807411287, -0.6010536646032245, -1.2662424075749568], [-0.3710451618282834, -0.11265730320380807, -2.315646431213898], [-0.24552532002048744, 1.8598848945507132, -1.4267659957399854], [1.577042443617165, 0.0, 0.0], [2.292718146893918, 1.3915527243580494, 0.0], [2.341079856722332, 2.059852692894922, 1.4165023767064753], [1.1453916410070462, 2.033172555891601, 1.9759511228307884], [3.202217992284648, 1.4422259961598178, 2.198403062836331], [2.721555559149519, 3.3152059779287892, 1.2745358845394053], [1.6292323391939894, 2.2122558673105726, -0.8090479336198849], [3.545586830094389, 1.26003922143106, -0.42809146886198096], [1.9974224573334833, -0.6906780683055249, 1.0535722235493006], [1.9277183224308936, -0.6529932317206335, -1.1102241252095322], [-0.3501493572534739, 0.7132786644586373, 1.070608697319931], [-0.42668432219275976, -1.2491488329668494, 0.15362238828850333], [-4.205197626772494, 4.07677799320728, 0.5117452221808406]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0408', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
