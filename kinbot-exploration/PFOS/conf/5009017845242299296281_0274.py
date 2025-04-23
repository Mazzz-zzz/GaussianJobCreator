import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0274'
logfile = 'conf/5009017845242299296281_0274.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863821, 0.6217394783082127, -1.2501828803165014], [-0.3466020415139028, 2.127818130564323, -1.5663863369811224], [1.162453549046707, 2.4914860914100156, -1.7945884035232595], [1.9243606525365986, 2.688842566201291, -0.438073329330714], [1.5424767514304083, 1.7679870357247445, 0.4295047485064743], [3.237845094852695, 2.6011438926710686, -0.627724596153851], [1.5819802655580881, 4.373096286345695, 0.3154042137368805], [1.9592615238581805, 4.3322211909891095, 1.685958512743339], [2.050177455734448, 5.339422692612445, -0.6078300140247241], [0.009894955985936404, 4.348157248361431, 0.22338001245885236], [1.7462981079693154, 1.5025680769231304, -2.4655267346281255], [1.253597860136186, 3.618508316135832, -2.4888098530759883], [-1.0117971004258963, 2.410996235613989, -2.6898422467352834], [-0.8215282967123811, 2.883031579863843, -0.5833314461381587], [-2.011954187959716, 0.6042709716797771, -1.0485616399675899], [-0.4099470658637771, -0.11587296658230463, -2.3181096973944277], [1.5770424436171684, 0.0, 0.0], [2.2927181468939164, 1.3915527243580577, 0.0], [1.6005215470082483, 2.4407219045638975, -0.9357086002340288], [1.376069276137164, 1.9145350254105704, -2.125714308102145], [0.46117091336274973, 2.8519467789832493, -0.41882897761344806], [2.4048649076934923, 3.47808001118291, -1.0689303403306873], [3.5400592233304646, 1.2295174124846084, -0.431521055692747], [2.308846803952251, 1.8960947387583764, 1.2310220414904687], [1.9974224573334873, -0.690678068305525, 1.0535722235492995], [1.9277183224308967, -0.652993231720626, -1.110224125209532], [-0.3501493572534735, 0.7132786644586341, 1.0706086973199345], [-0.42668432219275504, -1.2491488329668516, 0.15362238828850133], [-0.3710276935280501, 4.048170882852776, 1.0605589074939277]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0274', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
