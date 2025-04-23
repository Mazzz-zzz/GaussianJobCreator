import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0203'
logfile = 'conf/5009017845242299296281_0203.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586382, 0.6217394783082133, -1.250182880316501], [-0.34660204151390284, 2.127818130564324, -1.5663863369811217], [1.1624535490467065, 2.491486091410016, -1.79458840352326], [1.9243606525365986, 2.688842566201291, -0.438073329330714], [1.5424767514304085, 1.7679870357247447, 0.42950474850647435], [3.2378450948526947, 2.601143892671068, -0.6277245961538509], [1.5819802655580877, 4.373096286345695, 0.3154042137368825], [0.1921341302265273, 4.644568250446804, 0.18478376885551676], [2.288956239795058, 4.426539708800599, 1.5412449636289753], [2.36862059055278, 5.2530054928834, -0.7274281875392955], [1.746298107969314, 1.502568076923131, -2.4655267346281264], [1.2535978601361861, 3.618508316135833, -2.4888098530759866], [-1.0117971004258965, 2.4109962356139896, -2.6898422467352843], [-0.8215282967123819, 2.883031579863843, -0.5833314461381572], [-2.0119541879597156, 0.6042709716797787, -1.0485616399675928], [-0.40994706586377727, -0.11587296658230313, -2.3181096973944286], [1.5770424436171686, 0.0, 0.0], [2.2927181468939164, 1.3915527243580574, 0.0], [1.6005215470082481, 2.440721904563897, -0.9357086002340288], [1.3760692761371622, 1.9145350254105717, -2.125714308102145], [0.46117091336274996, 2.8519467789832493, -0.4188289776134476], [2.404864907693491, 3.478080011182911, -1.0689303403306871], [3.5400592233304646, 1.2295174124846087, -0.4315210556927487], [2.3088468039522514, 1.8960947387583766, 1.2310220414904691], [1.997422457333487, -0.690678068305526, 1.0535722235492986], [1.927718322430896, -0.6529932317206253, -1.1102241252095328], [-0.3501493572534738, 0.7132786644586349, 1.0706086973199334], [-0.42668432219275504, -1.2491488329668516, 0.15362238828850133], [1.7667412333238892, 5.586798730430628, -1.4073441386725756]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0203', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
