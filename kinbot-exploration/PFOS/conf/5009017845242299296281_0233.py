import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0233'
logfile = 'conf/5009017845242299296281_0233.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863807, 0.6217394783082059, -1.2501828803165083], [-2.270962283629192, 0.6501421835576421, -1.2334320314121858], [-2.9699917885127123, 1.3987418462127565, -0.04493060495240028], [-4.434530032795169, 1.8262065172130024, -0.4082305431339405], [-5.030583914311105, 0.8584013505576154, -1.0824522224071023], [-5.127712885957738, 2.088652764913048, 0.696088043009831], [-4.4665010348623895, 3.375387929270149, -1.4668514094871619], [-4.226968110733983, 4.495272964014938, -0.6236951398531391], [-3.721531756804504, 3.089411167353828, -2.636720991125125], [-5.992464116672815, 3.353489951035306, -1.8560811376810036], [-3.0171454311974584, 0.585971635842988, 1.0068168322315654], [-2.2839557356380094, 2.491572728979182, 0.26410822321819183], [-2.6212054717929205, 1.2630031390416123, -2.367697716836102], [-2.713377980741126, -0.6010536646032335, -1.266242407574963], [-0.37104516182827757, -0.11265730320381648, -2.315646431213901], [-0.24552532002048638, 1.8598848945507076, -1.4267659957399927], [1.577042443617164, 0.0, 0.0], [2.2927181468939155, 1.391552724358052, 0.0], [2.3410798567223248, 2.059852692894924, 1.4165023767064777], [1.145391641007037, 2.033172555891602, 1.975951122830787], [3.2022179922846385, 1.4422259961598192, 2.1984030628363356], [2.7215555591495106, 3.315205977928793, 1.274535884539405], [1.6292323391939887, 2.2122558673105743, -0.8090479336198868], [3.545586830094388, 1.2600392214310636, -0.4280914688619757], [1.9974224573334802, -0.69067806830552, 1.0535722235493066], [1.9277183224308976, -0.6529932317206335, -1.110224125209525], [-0.35014935725347807, 0.7132786644586384, 1.0706086973199265], [-0.42668432219275787, -1.2491488329668503, 0.1536223882885068], [-6.113344565466942, 3.676019135677995, -2.760141549848821]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0233', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
