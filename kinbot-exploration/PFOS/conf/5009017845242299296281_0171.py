import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0171'
logfile = 'conf/5009017845242299296281_0171.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863826, 0.7718203945763861, 1.1635336229088469], [-2.2709622836291934, 0.7431123812655647, 1.1797556627388974], [-3.0203184489305848, 1.3845311280592572, -0.04033723253932924], [-4.480058160568046, 0.8292457901518717, -0.1831774466796789], [-4.453761713276069, -0.3598415970703202, -0.7592691281731027], [-5.05670016854204, 0.73012625773328, 1.01120915467093], [-5.555811090753355, 1.9566623859997128, -1.228811600540189], [-6.700610364745076, 1.2211909591412446, -1.6422694706464422], [-5.612894691041406, 3.2047869210771136, -0.5622657437801356], [-4.587105203609392, 2.1176003861733994, -2.4601728367831184], [-3.0861165692439867, 2.7007512384555876, 0.13918095834253172], [-2.36392887439246, 1.1178871860085928, -1.1621587584935864], [-2.5961406643712652, -0.5518249635012433, 1.2269126406770385], [-2.6897770261666536, 1.3325843308485144, 2.2930463940297656], [-0.3710451618282835, 2.0617372872159114, 1.0602591291106065], [-0.24552532002049293, 0.3056731502912656, 2.324090564665856], [1.5770424436171668, 0.0, 0.0], [2.292718146893917, 1.3915527243580554, 0.0], [3.782335574419716, 1.3186147352454585, -0.48079377647244065], [4.419592497958775, 0.34925073248439464, 0.14976318527135185], [3.8473329977166415, 1.110443285289195, -1.7795740852228736], [4.36801041105158, 2.4684340853352067, -0.20560554420871807], [2.293896732220213, 1.8704189044736081, 1.2405689893126344], [1.649464944000881, 2.235214894314334, -0.8029305726284882], [1.9974224573334838, -0.6906780683055228, 1.0535722235493072], [1.9277183224308976, -0.6529932317206288, -1.110224125209524], [-0.3501493572534741, -1.2838136616209441, 0.0824130947386519], [-0.4266843221927539, 0.4915335256355426, -1.1586058166012276], [-4.816431986749273, 1.4829006644540734, -3.153374017163045]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0171', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
