import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0283'
logfile = 'conf/5009017845242299296281_0283.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863818, 0.7718203945763839, 1.1635336229088502], [-2.2709622836291907, 0.7431123812655646, 1.179755662738902], [-3.0203184489305834, 1.3845311280592618, -0.040337232539322], [-4.480058160568046, 0.8292457901518818, -0.18317744667967145], [-5.202953426281215, 1.6517668770716232, -0.9230311075200041], [-4.464154537777424, -0.37536056335173595, -0.7467115170195916], [-5.325329553432012, 0.6531877726237025, 1.4830107713751979], [-6.725522199345814, 0.5396874740532336, 1.2616304552598758], [-4.563283221999027, -0.27393512404454634, 2.234732325456914], [-5.018129818548461, 2.0912608638679253, 2.0470291548838886], [-3.0861165692439863, 2.700751238455591, 0.13918095834254274], [-2.3639288743924602, 1.1178871860085973, -1.16215875849358], [-2.5961406643712626, -0.5518249635012418, 1.2269126406770408], [-2.6897770261666487, 1.332584330848513, 2.2930463940297723], [-0.3710451618282798, 2.0617372872159123, 1.0602591291106105], [-0.24552532002048877, 0.30567315029126363, 2.3240905646658563], [1.577042443617165, 0.0, 0.0], [2.2927181468939173, 1.3915527243580539, 0.0], [1.6005215470082566, 2.4407219045638966, -0.9357086002340287], [1.376069276137168, 1.9145350254105806, -2.1257143081021423], [0.4611709133627595, 2.8519467789832564, -0.4188289776134503], [2.404864907693508, 3.4780800111829087, -1.068930340330682], [3.540059223330469, 1.2295174124846036, -0.43152105569274357], [2.3088468039522505, 1.89609473875837, 1.2310220414904727], [1.9974224573334856, -0.6906780683055236, 1.0535722235493008], [1.9277183224308962, -0.6529932317206264, -1.1102241252095315], [-0.35014935725347524, -1.2838136616209428, 0.08241309473865066], [-0.4266843221927563, 0.49153352563554625, -1.1586058166012247], [-4.233362975845368, 2.4599766596281363, 1.6178997440452285]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0283', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
