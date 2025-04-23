import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0477'
logfile = 'conf/5009017845242299296281_0477.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863781, 0.7718203945763846, 1.1635336229088538], [-2.2709622836291867, 0.7431123812655684, 1.1797556627389083], [-2.9699917885127096, -0.660459877810206, 1.2338112746327947], [-4.434530032795165, -0.5595652376518028, 1.7856565080301383], [-5.118760473256557, -1.6369672213233482, 1.4428218350409119], [-4.423212373540424, -0.4450548041217596, 3.1106694557988646], [-5.334436835556407, 0.9399596749116969, 1.1050273034991618], [-5.0431109608800355, 1.0412804304034988, -0.2832429681196695], [-6.6450608998698, 0.9110536536080365, 1.6405511364170382], [-4.530265242807154, 2.055641721998861, 1.8725585725895624], [-3.0171454311974584, -1.1649147715918018, 0.0040579064213636575], [-2.2839557356380067, -1.4745107951449332, 2.025711167063388], [-2.621205471792916, 1.4189848017416533, 2.2776416618875945], [-2.7133779807411256, 1.3971249246107036, 0.11259346120336604], [-0.3710451618282739, 2.0617372872159114, 1.0602591291106191], [-0.24552532002048283, 0.3056731502912572, 2.3240905646658554], [1.5770424436171655, 0.0, 0.0], [2.2927181468939204, 1.3915527243580528, 0.0], [1.6005215470082534, 2.4407219045638975, -0.9357086002340258], [1.3760692761371636, 1.9145350254105828, -2.1257143081021406], [0.46117091336276417, 2.851946778983258, -0.4188289776134443], [2.4048649076935082, 3.4780800111829073, -1.0689303403306845], [3.540059223330463, 1.2295174124845976, -0.43152105569274823], [2.308846803952263, 1.8960947387583666, 1.231022041490469], [1.9974224573334833, -0.6906780683055344, 1.0535722235492946], [1.9277183224308883, -0.6529932317206242, -1.1102241252095428], [-0.3501493572534812, -1.2838136616209432, 0.08241309473864497], [-0.42668432219275787, 0.4915335256355536, -1.1586058166012236], [-4.141293426963926, 1.6935296680348217, 2.6809756515466896]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0477', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
