import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0257'
logfile = 'conf/5009017845242299296281_0257.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863898, 0.6217394783082143, -1.2501828803165], [-0.3976197158559613, -0.0756648590159576, -2.633910119820667], [-1.1233882121466392, 0.4940238789133657, -3.902883790150537], [-0.3678111583333681, 0.1149613201546942, -5.223868290930236], [-1.1821704538825721, 0.23639236066085137, -6.257587490341371], [0.6879514164795726, 0.9051956426165999, -5.3964809839181385], [0.26055307355107504, -1.6528948974501971, -5.185059566666355], [1.4264341127623, -1.692774024579163, -4.371694909555421], [-0.8683369306174475, -2.486431734998271, -4.994847919844323], [0.6934656418909404, -1.773642166979223, -6.694546047190003], [-2.35153528541847, -0.013330338437769183, -3.959787366978751], [-1.1885301103712445, 1.8171105205812574, -3.827863304783672], [0.918670262921738, 0.05926285465343761, -2.8183083047232302], [-0.6863393184328456, -1.366668279436513, -2.5232286857846633], [-0.25604457595343505, 1.8735740976390427, -1.3916701657561914], [-2.0076024771874525, 0.6427130616946727, -1.0543092166280603], [1.5770424436171635, 0.0, 0.0], [2.292718146893908, 1.3915527243580599, 0.0], [1.6005215470082343, 2.4407219045639, -0.9357086002340294], [1.376069276137151, 1.9145350254105722, -2.125714308102147], [0.4611709133627384, 2.851946778983245, -0.4188289776134513], [2.40486490769348, 3.4780800111829144, -1.0689303403306867], [3.5400592233304584, 1.2295174124846227, -0.431521055692748], [2.308846803952242, 1.8960947387583804, 1.231022041490467], [1.997422457333483, -0.6906780683055216, 1.0535722235492988], [1.9277183224308934, -0.652993231720622, -1.1102241252095313], [-0.35014935725348284, 0.7132786644586324, 1.070608697319934], [-0.42668432219275554, -1.249148832966853, 0.15362238828850264], [0.1970420147741041, -1.1456221188817794, -7.2377925734778845]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0257', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
