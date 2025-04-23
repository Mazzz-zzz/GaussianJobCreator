import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0014'
logfile = 'conf/5009017845242299296281_0014.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863861, -1.393559872884596, 0.08664925740765192], [-2.2709622836291956, -1.3932545648232129, 0.053676368673271066], [-2.9970239643019054, -2.7837283183799193, 0.08524408716293146], [-3.0516123514917433, -3.4433745877869066, -1.3365407299748304], [-4.021172178383429, -2.896073648552945, -2.048522026827762], [-1.892287467807033, -3.2800720040577174, -1.967570143528445], [-3.3760113189425556, -5.289049349093572, -1.2369750907842632], [-3.8039531646981053, -5.7351505791457615, -2.517713345466792], [-2.2964465318638347, -5.858101018474394, -0.518590208062902], [-4.621512289586947, -5.249481741506512, -0.2738009806864454], [-2.3329497198238265, -3.5970575181208235, 0.901634044950072], [-4.241488402883609, -2.6316641353094163, 0.5195418500188193], [-2.6427825570536205, -0.7109480038655779, 1.1404306775613242], [-2.667789377892034, -0.7313203152515839, -1.0264488679511474], [-0.37104516182828834, -1.9490799840121078, 1.2553873021032877], [-0.2455253200204925, -2.165558044841978, -0.8973245689258773], [1.577042443617164, 0.0, 0.0], [2.2927181468939164, 1.3915527243580539, 0.0], [3.782335574419716, 1.318614735245457, -0.48079377647244037], [4.419592497958777, 0.34925073248439087, 0.14976318527136112], [3.847332997716635, 1.1104432852891999, -1.77957408522288], [4.368010411051583, 2.4684340853352014, -0.2056055442087094], [2.293896732220213, 1.8704189044736084, 1.2405689893126326], [1.6494649440008786, 2.235214894314339, -0.8029305726284883], [1.9974224573334787, -0.6906780683055256, 1.053572223549302], [1.9277183224308934, -0.6529932317206265, -1.1102241252095324], [-0.35014935725347257, 0.5705349971623099, -1.1530217920585821], [-0.42668432219276103, 0.757615307331307, 1.004983428312724], [-5.442540181202505, -5.2714394378063965, -0.7850674938266146]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0014', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
