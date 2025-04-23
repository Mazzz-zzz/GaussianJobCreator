import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0344'
logfile = 'conf/5009017845242299296281_0344.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863809, 0.6217394783082097, -1.2501828803165052], [-0.3466020415138996, 2.12781813056432, -1.566386336981132], [-0.7363023803695162, 3.2026189464980015, -0.4919021834050506], [-2.245266873553545, 3.61521560384561, -0.603146579477167], [-2.9870600258534608, 2.5464156154084763, -0.835334966341661], [-2.6496830196655057, 4.199339702913139, 0.5211961877453115], [-2.5333983411176324, 4.835652692402317, -1.9992758361685343], [-2.09647104819226, 6.11789670754432, -1.5664065324065186], [-2.094654062919407, 4.208503642557027, -3.190722141830052], [-4.108201473709184, 4.81248651985339, -2.000301427835373], [-0.5266683954723083, 2.6951529346344745, 0.7195116152968816], [0.004125934750531798, 4.291386672281076, -0.6555498260725718], [0.9809271994195118, 2.1566775839163728, -1.7135621949231852], [-0.9145558447621632, 2.455752201920101, -2.720590366825878], [-2.0119541879597147, 0.6042709716797786, -1.0485616399676003], [-0.40994706586377255, -0.11587296658230949, -2.3181096973944317], [1.5770424436171668, 0.0, 0.0], [2.292718146893915, 1.391552724358055, 0.0], [1.6005215470082443, 2.4407219045638917, -0.9357086002340341], [1.376069276137165, 1.9145350254105638, -2.1257143081021477], [0.4611709133627504, 2.851946778983248, -0.4188289776134555], [2.40486490769349, 3.478080011182903, -1.0689303403306918], [3.5400592233304677, 1.2295174124846078, -0.43152105569274435], [2.308846803952245, 1.8960947387583778, 1.231022041490468], [1.997422457333481, -0.6906780683055256, 1.0535722235493037], [1.927718322430896, -0.6529932317206297, -1.1102241252095266], [-0.3501493572534781, 0.7132786644586365, 1.0706086973199278], [-0.42668432219275715, -1.2491488329668499, 0.15362238828849772], [-4.455412368011851, 5.5447313633935424, -1.4718844800730977]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0344', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
