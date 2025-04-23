import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0385'
logfile = 'conf/5009017845242299296281_0385.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863796, 0.7718203945763868, 1.1635336229088502], [-2.2709622836291894, 0.7431123812655711, 1.1797556627389032], [-2.969991788512711, -0.6604598778102055, 1.233811274632794], [-4.434530032795166, -0.5595652376518016, 1.7856565080301345], [-5.118760473256556, -1.636967221323349, 1.44282183504091], [-4.423212373540425, -0.44505480412175386, 3.1106694557988632], [-5.334436835556408, 0.9399596749116975, 1.1050273034991547], [-6.728895866558153, 0.754355760193769, 1.3135128078131466], [-4.613374145520476, 2.0813000197186304, 1.5325405617756103], [-5.018585672357864, 0.7209459420014732, -0.4223278147729051], [-3.017145431197457, -1.1649147715918058, 0.004057906421363671], [-2.2839557356380076, -1.4745107951449312, 2.025711167063389], [-2.621205471792919, 1.4189848017416593, 2.277641661887586], [-2.713377980741125, 1.3971249246107038, 0.11259346120335982], [-0.3710451618282775, 2.0617372872159136, 1.0602591291106136], [-0.24552532002048547, 0.30567315029126335, 2.3240905646658545], [1.5770424436171655, 0.0, 0.0], [2.29271814689392, 1.3915527243580548, 0.0], [1.6005215470082519, 2.440721904563897, -0.9357086002340271], [1.376069276137166, 1.9145350254105806, -2.1257143081021423], [0.4611709133627555, 2.8519467789832547, -0.4188289776134464], [2.4048649076935, 3.4780800111829078, -1.0689303403306831], [3.5400592233304637, 1.2295174124846027, -0.43152105569274535], [2.3088468039522545, 1.8960947387583689, 1.2310220414904698], [1.9974224573334836, -0.6906780683055318, 1.0535722235492957], [1.9277183224308887, -0.6529932317206257, -1.1102241252095406], [-0.3501493572534805, -1.2838136616209432, 0.0824130947386484], [-0.42668432219275826, 0.4915335256355488, -1.158605816601228], [-4.889526728848133, 1.5718645749678157, -0.8641916579042966]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0385', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
