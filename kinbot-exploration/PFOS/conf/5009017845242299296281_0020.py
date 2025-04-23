import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0020'
logfile = 'conf/5009017845242299296281_0020.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863819, 0.6217394783082155, -1.2501828803165027], [-2.270962283629193, 0.650142183557653, -1.233432031412176], [-3.020318448930584, -0.7271986321270438, -1.1788705129599821], [-4.480058160568045, -0.573259217300906, -0.6265591969129879], [-4.45376171327607, -0.47762555477201035, 0.6912665284878049], [-5.056700168542037, 0.5106696876177758, -1.137912464502551], [-5.555811090753356, -2.0425132555326706, -1.0801135326351394], [-5.935868144708009, -1.9112395063211098, -2.4442084143487435], [-4.931749541331282, -3.192089888270275, -0.5375478503029378], [-6.7987732354635675, -1.7025466603097759, -0.17456924837349033], [-3.086116569243985, -1.2298413735800922, -2.408509660976095], [-2.3639288743924594, -1.5654026010903246, -0.38703932240175537], [-2.596140664371265, 1.3384499968011858, -0.13556188350402096], [-2.6897770261666514, 1.3195442638618284, -2.3005750802147826], [-0.3710451618282821, -0.11265730320379824, -2.3156464312138976], [-0.2455253200204908, 1.8598848945507216, -1.4267659957399774], [1.5770424436171673, 0.0, 0.0], [2.292718146893916, 1.3915527243580605, 0.0], [1.6005215470082446, 2.4407219045638984, -0.9357086002340301], [1.3760692761371636, 1.9145350254105726, -2.125714308102145], [0.46117091336274996, 2.8519467789832516, -0.41882897761345195], [2.4048649076934945, 3.4780800111829127, -1.0689303403306878], [3.5400592233304646, 1.2295174124846155, -0.4315210556927484], [2.308846803952254, 1.8960947387583786, 1.2310220414904691], [1.997422457333487, -0.6906780683055237, 1.053572223549291], [1.9277183224308967, -0.652993231720623, -1.1102241252095322], [-0.3501493572534739, 0.713278664458636, 1.0706086973199331], [-0.42668432219275554, -1.2491488329668492, 0.15362238828850217], [-7.468324495646526, -1.2282948166259808, -0.6871594213635064]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0020', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
