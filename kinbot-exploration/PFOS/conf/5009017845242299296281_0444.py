import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0444'
logfile = 'conf/5009017845242299296281_0444.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863843, 0.6217394783082139, -1.250182880316502], [-0.3976197158559566, -0.07566485901595763, -2.633910119820668], [1.0879414097563251, -0.08367235852192167, -3.138512306045848], [2.100632021959493, -0.3015263760374457, -1.960979924321287], [3.2636602341713647, -0.708743058804942, -2.438517246372761], [2.2719309703990027, 0.8294205503577342, -1.2823668088469176], [1.48866165575111, -1.5931392216697093, -0.7448829466934757], [2.5928813243933715, -2.0212813565520817, 0.04240622732833154], [0.2718090320333246, -1.107780636329468, -0.20724598048332032], [1.134464322195662, -2.728240045041822, -1.7776622065022532], [1.2395907752722295, -1.0741885454775717, -4.013081746298286], [1.3717373297746493, 1.0745873687029925, -3.72016863203124], [-0.7843720924895169, -1.3445493546988805, -2.4752210448161542], [-1.1581983763610537, 0.49302227342429195, -3.56149530020636], [-0.25604457595343105, 1.8735740976390431, -1.3916701657561947], [-2.007602477187449, 0.642713061694673, -1.054309216628062], [1.5770424436171655, 0.0, 0.0], [2.2927181468939124, 1.3915527243580588, 0.0], [1.6005215470082415, 2.4407219045638966, -0.935708600234033], [1.376069276137155, 1.9145350254105709, -2.125714308102144], [0.4611709133627433, 2.8519467789832476, -0.41882897761344884], [2.4048649076934887, 3.4780800111829118, -1.068930340330686], [3.540059223330462, 1.2295174124846147, -0.4315210556927462], [2.3088468039522456, 1.8960947387583778, 1.231022041490468], [1.997422457333486, -0.6906780683055225, 1.0535722235492995], [1.9277183224308934, -0.6529932317206256, -1.110224125209532], [-0.350149357253477, 0.7132786644586341, 1.0706086973199302], [-0.42668432219275615, -1.2491488329668503, 0.1536223882885023], [0.9564676515766036, -2.34653264220792, -2.648626632048057]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0444', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
