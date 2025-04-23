import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0215'
logfile = 'conf/5009017845242299296281_0215.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863828, 0.6217394783082106, -1.2501828803165038], [-2.270962283629193, 0.6501421835576464, -1.2334320314121807], [-2.969991788512713, 1.3987418462127557, -0.04493060495238967], [-2.1956407956954176, 2.701984397794214, 0.35697056532746907], [-2.9912731058630335, 3.5033560351048907, 1.0434650438255453], [-1.1300133907845045, 2.3990794081170193, 1.0928964051269063], [-1.5840214570309339, 3.6417954560905574, -1.1478019392028767], [-1.2617908272916094, 4.966762823772167, -0.7439951739714586], [-0.6920488892418214, 2.7826859433857702, -1.834586152140139], [-2.933517072785173, 3.675195743736788, -1.9591518228253695], [-4.2015349880113835, 1.7379422988872626, -0.41538431755355903], [-3.0244357979279366, 0.6050007412875016, 1.0168795050389554], [-2.621205471792927, 1.2630031390416194, -2.367697716836093], [-2.713377980741127, -0.6010536646032306, -1.2662424075749596], [-0.3710451618282802, -0.11265730320380615, -2.315646431213901], [-0.24552532002049016, 1.8598848945507165, -1.4267659957399828], [1.5770424436171655, 0.0, 0.0], [2.2927181468939137, 1.3915527243580577, 0.0], [1.600521547008244, 2.4407219045638984, -0.9357086002340305], [1.376069276137164, 1.9145350254105729, -2.1257143081021463], [0.46117091336274796, 2.8519467789832498, -0.41882897761345195], [2.4048649076934954, 3.4780800111829135, -1.0689303403306878], [3.540059223330464, 1.2295174124846138, -0.43152105569274324], [2.308846803952247, 1.8960947387583769, 1.2310220414904693], [1.9974224573334836, -0.6906780683055243, 1.0535722235493021], [1.9277183224308978, -0.6529932317206243, -1.1102241252095282], [-0.3501493572534782, 0.7132786644586352, 1.070608697319932], [-0.4266843221927564, -1.2491488329668519, 0.1536223882885025], [-2.75764586516072, 3.5968192019828233, -2.907249298155522]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0215', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
