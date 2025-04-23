import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0238'
logfile = 'conf/5009017845242299296281_0238.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863843, 0.6217394783082144, -1.2501828803165003], [-2.2709622836291947, 0.6501421835576481, -1.2334320314121732], [-2.9699917885127136, 1.398741846212756, -0.04493060495238037], [-3.056739804078713, 0.500064566203076, 1.2374262436633878], [-1.9299416078656424, -0.17317330177087917, 1.3909793003371522], [-3.2751648795012236, 1.2495736306950789, 2.314193019350682], [-4.462079523721017, -0.7386873768807439, 1.127075770095076], [-4.239230824117429, -1.7519461414208217, 2.0997402245590506], [-5.660645840812833, 0.007769908164212347, 1.0193860454956845], [-4.136425722847151, -1.3258498871817905, -0.2976110089136924], [-2.2659422418234905, 2.4877087631800117, 0.2508004249777433], [-4.202897212097874, 1.7436062796993266, -0.3934166742568534], [-2.6212054717929307, 1.2630031390416228, -2.3676977168360835], [-2.71337798074113, -0.6010536646032287, -1.266242407574954], [-0.37104516182828434, -0.11265730320380446, -2.315646431213898], [-0.2455253200204931, 1.8598848945507178, -1.4267659957399785], [1.5770424436171637, 0.0, 0.0], [2.2927181468939146, 1.3915527243580585, 0.0], [3.7823355744197134, 1.318614735245458, -0.4807937764724519], [4.419592497958775, 0.3492507324843943, 0.14976318527135069], [3.847332997716635, 1.1104432852892012, -1.7795740852228812], [4.368010411051579, 2.4684340853352102, -0.20560554420872212], [2.2938967322202166, 1.8704189044736077, 1.2405689893126286], [1.649464944000877, 2.2352148943143373, -0.8029305726284892], [1.997422457333485, -0.6906780683055207, 1.0535722235493017], [1.9277183224308962, -0.6529932317206228, -1.1102241252095313], [-0.35014935725347474, 0.7132786644586353, 1.0706086973199362], [-0.4266843221927538, -1.2491488329668503, 0.15362238828850455], [-3.6343716860935866, -0.6838820826267326, -0.8189376519986283]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0238', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
