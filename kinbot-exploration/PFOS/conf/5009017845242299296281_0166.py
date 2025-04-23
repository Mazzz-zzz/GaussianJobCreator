import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0166'
logfile = 'conf/5009017845242299296281_0166.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863829, 0.7718203945763823, 1.1635336229088502], [-2.2709622836291947, 0.7431123812655587, 1.179755662738902], [-2.997023964301903, 1.4656877041954754, 2.3681573973696826], [-3.0516123514917357, 0.5642090685426548, 3.650320232756651], [-3.2793160186248658, 1.3122077232124991, 4.715754883809086], [-4.013520036470854, -0.3470252206198689, 3.535056550767916], [-1.4404858172046033, -0.35858040705153127, 3.9230732444055194], [-1.3958302980578645, -1.4531530551280987, 3.0162394059861595], [-0.4216488399427635, 0.6181634544205777, 4.038396331810095], [-1.736726766178271, -0.9114086423563111, 5.367775591676397], [-2.332949719823824, 2.579366746904097, 2.6643261670914025], [-4.241488402883606, 1.7657685081001735, 2.019317070396956], [-2.642782557053619, 1.3431159399559982, 0.04548369333676725], [-2.6677893778920323, -0.5232706377056777, 1.1465664052870914], [-0.37104516182828845, 2.0617372872159105, 1.0602591291106107], [-0.24552532002048855, 0.30567315029126124, 2.3240905646658545], [1.5770424436171655, 0.0, 0.0], [2.2927181468939133, 1.3915527243580603, 0.0], [2.3410798567223217, 2.0598526928949306, 1.4165023767064722], [1.145391641007034, 2.0331725558915914, 1.9759511228307955], [3.202217992284642, 1.4422259961598332, 2.198403062836327], [2.7215555591495013, 3.315205977928796, 1.2745358845394066], [1.6292323391939745, 2.2122558673105797, -0.8090479336198859], [3.54558683009438, 1.2600392214310843, -0.42809146886198324], [1.9974224573334856, -0.6906780683055223, 1.053572223549295], [1.927718322430893, -0.6529932317206273, -1.1102241252095362], [-0.35014935725347374, -1.2838136616209457, 0.08241309473864741], [-0.4266843221927599, 0.4915335256355461, -1.1586058166012274], [-0.9262218000357515, -0.9230104746980652, 5.895894407891676]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0166', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
