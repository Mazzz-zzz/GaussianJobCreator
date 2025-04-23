import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0267'
logfile = 'conf/5009017845242299296281_0267.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863781, -1.3935598728846033, 0.08664925740764866], [-2.270962283629187, -1.3932545648232253, 0.05367636867327403], [-2.9699917885127096, -0.7382819684025533, -1.1888806696804035], [-4.434530032795165, -1.2666412795612068, -1.3774259648962097], [-5.118760473256557, -0.43103675161863103, -2.1390661163489], [-4.423212373540424, -2.471391369437247, -1.9407634943451886], [-5.334436835556407, -1.4269615541615401, 0.2615153052569081], [-5.043110960880036, -0.2753446093668123, 1.0433967892528582], [-6.6450608998698, -1.8762857871486034, -0.03127995997333542], [-4.530265242807154, -2.649504154936319, 0.8439586660354184], [-3.017145431197458, 0.5789431357488213, -1.0108747386529282], [-2.2839557356380067, -1.0170619338342468, -2.2898193902816004], [-2.621205471792916, -2.6819879407832907, 0.09005605494849686], [-2.7133779807411256, -0.7960712600074867, 1.1536489463716033], [-0.3710451618282739, -1.949079984012119, 1.25538730210328], [-0.24552532002048283, -2.165558044841979, -0.8973245689258833], [1.5770424436171655, 0.0, 0.0], [2.292718146893915, 1.3915527243580563, 0.0], [3.78233557441971, 1.318614735245477, -0.48079377647243515], [4.419592497958773, 0.34925073248440774, 0.14976318527134913], [3.8473329977166273, 1.1104432852892214, -1.7795740852228716], [4.368010411051572, 2.46843408533521, -0.2056055442087133], [2.2938967322202135, 1.870418904473608, 1.2405689893126364], [1.6494649440008702, 2.235214894314346, -0.8029305726284743], [1.9974224573334871, -0.6906780683055235, 1.0535722235492941], [1.9277183224308956, -0.6529932317206135, -1.110224125209536], [-0.35014935725347907, 0.5705349971623113, -1.153021792058573], [-0.42668432219275787, 0.7576153073312951, 1.0049834283127324], [-4.1412934269639265, -3.1685578551843814, 0.12615188880743172]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0267', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
