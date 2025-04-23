import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0000'
logfile = 'conf/5009017845242299296281_0000.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863796, 0.7718203945763865, 1.163533622908852], [-2.2709622836291894, 0.7431123812655698, 1.1797556627389056], [-2.9699917885127096, -0.6604598778102079, 1.2338112746327956], [-4.434530032795166, -0.5595652376518061, 1.7856565080301376], [-4.4121494903032055, -0.46329971437302003, 3.1035085495820818], [-5.049857983877611, 0.5028813412615518, 1.2743065103524205], [-5.456246996148598, -2.0667053949410272, 1.3314689650514], [-5.840354017548495, -1.9490077204655891, -0.032730456347233254], [-4.791164773781588, -3.193046400996997, 1.874035410029529], [-6.711000276427691, -1.7721280585589976, 2.236654597023056], [-3.0171454311974557, -1.1649147715918085, 0.004057906421367819], [-2.2839557356380054, -1.4745107951449323, 2.0257111670633905], [-2.621205471792919, 1.4189848017416555, 2.2776416618875888], [-2.713377980741127, 1.397124924610699, 0.11259346120336129], [-0.3710451618282788, 2.0617372872159145, 1.0602591291106138], [-0.24552532002048513, 0.3056731502912657, 2.3240905646658563], [1.5770424436171648, 0.0, 0.0], [2.292718146893918, 1.3915527243580554, 0.0], [3.782335574419717, 1.3186147352454523, -0.4807937764724478], [4.41959249795878, 0.3492507324843904, 0.1497631852713477], [3.847332997716636, 1.110443285289191, -1.7795740852228836], [4.368010411051585, 2.4684340853351996, -0.20560554420872246], [2.2938967322202215, 1.8704189044736006, 1.2405689893126302], [1.6494649440008804, 2.2352148943143373, -0.8029305726284887], [1.9974224573334836, -0.6906780683055285, 1.053572223549297], [1.9277183224308905, -0.6529932317206245, -1.1102241252095384], [-0.35014935725347873, -1.2838136616209455, 0.08241309473865083], [-0.42668432219275915, 0.4915335256355483, -1.158605816601224], [-7.397108553229822, -1.3223216833423372, 1.723924605688382]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0000', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
