import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0417'
logfile = 'conf/5009017845242299296281_0417.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863828, 0.7718203945763874, 1.1635336229088442], [-0.39761971585595585, 2.3188655045575945, 1.25142736982874], [1.0879414097563267, 2.759867566386749, 1.496793764948375], [1.473042089945535, 2.685713952893136, 3.0152745408311827], [0.9614495666061577, 1.593462181345954, 3.5554236875418703], [2.7953494776680783, 2.6726461892025375, 3.1575153453829032], [0.8248630198382828, 4.164243062182728, 3.9720811551337913], [1.6664714540296413, 5.276724521658703, 3.695474478370246], [-0.5819280292581736, 4.174128097515165, 3.810223186820319], [1.1337571966355793, 3.6467518733462767, 5.4271859234671105], [1.894753649284141, 1.9494330700686118, 0.8176461448894697], [1.2598103705638592, 4.009816458287683, 1.086283558698185], [-0.7843720924895148, 2.815878982142089, 0.07319662459686575], [-1.1581983763610528, 2.837834268725452, 2.207717463520164], [-0.25604457595342933, 0.268434668414249, 2.318397847306013], [-2.007602477187446, 0.591702034196638, 1.0837604470856892], [1.5770424436171653, 0.0, 0.0], [2.2927181468939186, 1.3915527243580543, 0.0], [1.6005215470082568, 2.440721904563902, -0.9357086002340249], [1.3760692761371742, 1.9145350254105775, -2.125714308102145], [0.4611709133627624, 2.851946778983257, -0.41882897761344895], [2.404864907693505, 3.478080011182911, -1.0689303403306785], [3.5400592233304673, 1.2295174124846053, -0.43152105569274], [2.3088468039522487, 1.8960947387583729, 1.2310220414904722], [1.9974224573334847, -0.6906780683055254, 1.0535722235493008], [1.9277183224308967, -0.6529932317206304, -1.1102241252095286], [-0.35014935725347807, -1.2838136616209428, 0.08241309473865181], [-0.4266843221927558, 0.4915335256355436, -1.1586058166012314], [1.9831461677131388, 3.99368839639882, 5.734001521962399]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0417', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
