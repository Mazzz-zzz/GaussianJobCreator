import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0459'
logfile = 'conf/5009017845242299296281_0459.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863847, 0.7718203945763806, 1.1635336229088478], [-2.270962283629195, 0.7431123812655602, 1.1797556627388996], [-2.9970239643019054, 1.4656877041954737, 2.3681573973696795], [-3.0516123514917393, 0.5642090685426543, 3.6503202327566475], [-3.279316018624869, 1.3122077232125022, 4.715754883809081], [-4.013520036470855, -0.34702522061986874, 3.5350565507679144], [-1.440485817204605, -0.3585804070515311, 3.9230732444055176], [-1.4139902731590983, -0.819217681591184, 5.268265533730883], [-1.2483498686059542, -1.1927273355762995, 2.79496001846727], [-0.4537598945252275, 0.8614753930789404, 3.7875061839168978], [-2.3329497198238234, 2.579366746904097, 2.664326167091402], [-4.241488402883606, 1.7657685081001753, 2.0193170703969545], [-2.642782557053619, 1.3431159399559982, 0.045483693336765454], [-2.667789377892034, -0.5232706377056795, 1.1465664052870903], [-0.371045161828288, 2.061737287215911, 1.060259129110611], [-0.24552532002049055, 0.3056731502912611, 2.324090564665853], [1.5770424436171655, 0.0, 0.0], [2.2927181468939137, 1.3915527243580585, 0.0], [3.782335574419707, 1.3186147352454658, -0.48079377647244814], [4.419592497958774, 0.3492507324843981, 0.1497631852713489], [3.8473329977166357, 1.1104432852892028, -1.7795740852228747], [4.3680104110515705, 2.4684340853352125, -0.20560554420872612], [2.293896732220212, 1.8704189044736117, 1.2405689893126304], [1.6494649440008762, 2.2352148943143417, -0.8029305726284882], [1.9974224573334842, -0.6906780683055229, 1.0535722235492966], [1.9277183224308945, -0.6529932317206271, -1.1102241252095344], [-0.3501493572534743, -1.2838136616209452, 0.08241309473864625], [-0.4266843221927599, 0.4915335256355442, -1.1586058166012294], [-0.2384142873976419, 1.2181519118401458, 4.660644498934191]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0459', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
