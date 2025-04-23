import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0083'
logfile = 'conf/5009017845242299296281_0083.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863796, 0.7718203945763865, 1.163533622908852], [-2.2709622836291894, 0.7431123812655698, 1.1797556627389056], [-2.9699917885127096, -0.6604598778102079, 1.2338112746327956], [-4.434530032795166, -0.5595652376518061, 1.7856565080301376], [-4.4121494903032055, -0.46329971437302003, 3.1035085495820818], [-5.049857983877611, 0.5028813412615518, 1.2743065103524205], [-5.456246996148598, -2.0667053949410272, 1.3314689650514], [-5.840354017548495, -1.9490077204655891, -0.032730456347233254], [-4.791164773781588, -3.193046400996997, 1.874035410029529], [-6.711000276427691, -1.7721280585589976, 2.236654597023056], [-3.0171454311974557, -1.1649147715918085, 0.004057906421367819], [-2.2839557356380054, -1.4745107951449323, 2.0257111670633905], [-2.621205471792919, 1.4189848017416555, 2.2776416618875888], [-2.713377980741127, 1.397124924610699, 0.11259346120336129], [-0.3710451618282788, 2.0617372872159145, 1.0602591291106138], [-0.24552532002048513, 0.3056731502912657, 2.3240905646658563], [1.5770424436171648, 0.0, 0.0], [2.292718146893918, 1.3915527243580554, 0.0], [2.3410798567223257, 2.0598526928949235, 1.4165023767064777], [1.1453916410070435, 2.033172555891597, 1.9759511228307978], [3.202217992284652, 1.4422259961598183, 2.1984030628363294], [2.721555559149509, 3.3152059779287932, 1.2745358845394077], [1.6292323391939774, 2.2122558673105743, -0.8090479336198815], [3.5455868300943836, 1.2600392214310692, -0.4280914688619849], [1.9974224573334833, -0.6906780683055285, 1.0535722235492972], [1.9277183224308907, -0.6529932317206247, -1.1102241252095388], [-0.35014935725347873, -1.2838136616209455, 0.08241309473865083], [-0.42668432219275915, 0.4915335256355483, -1.158605816601224], [-7.076477711140584, -2.5980024687152214, 2.5835264610604782]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0083', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
