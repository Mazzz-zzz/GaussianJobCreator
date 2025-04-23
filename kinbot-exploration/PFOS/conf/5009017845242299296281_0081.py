import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0081'
logfile = 'conf/5009017845242299296281_0081.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863868, 0.7718203945763812, 1.1635336229088509], [-0.34660204151390495, 0.2926212946843381, 2.625937724192381], [1.1624535490467045, 0.3084161010830952, 3.0549844500983028], [1.3229555267971613, 0.3508819546540507, 4.6144782543247524], [0.43579710458927745, -0.4495079334091554, 5.179126732644479], [2.548490142211768, -0.029320092044959657, 4.964351589489722], [1.0607552256216175, 2.0797251118369364, 5.2956652118043905], [2.2341818168842673, 2.841309499474007, 5.039754557014124], [-0.24309140467177714, 2.4788600880503737, 4.913523867425547], [1.0184448819506249, 1.7166278846626593, 6.827628941220265], [1.7500514572469164, -0.7950580504033868, 2.6010405683398545], [1.763328214900581, 1.376908489527516, 2.5473394429719396], [-1.0117971004259, 1.1239736000383675, 3.4329051118380107], [-0.8215282967123799, -0.93633593874997, 2.788444311143953], [-2.011954187959719, 0.6059455318059125, 1.0475948322279962], [-0.4099470658637843, 2.065478369993783, 1.0587059160250796], [1.5770424436171646, 0.0, 0.0], [2.292718146893908, 1.3915527243580594, 0.0], [1.6005215470082366, 2.4407219045638993, -0.9357086002340314], [1.3760692761371505, 1.9145350254105713, -2.1257143081021423], [0.46117091336273663, 2.851946778983246, -0.41882897761345383], [2.40486490769348, 3.4780800111829153, -1.0689303403306913], [3.54005922333046, 1.2295174124846202, -0.43152105569275023], [2.308846803952243, 1.8960947387583824, 1.2310220414904627], [1.9974224573334844, -0.6906780683055244, 1.0535722235492986], [1.9277183224308965, -0.6529932317206228, -1.1102241252095326], [-0.35014935725347496, -1.2838136616209441, 0.08241309473864847], [-0.42668432219276103, 0.4915335256355452, -1.1586058166012252], [0.7519969798805651, 0.7943045879993396, 6.947159851911034]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0081', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
