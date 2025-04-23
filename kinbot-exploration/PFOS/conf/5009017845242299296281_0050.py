import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0050'
logfile = 'conf/5009017845242299296281_0050.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863873, -1.3935598728845933, 0.08664925740765798], [-2.2709622836291974, -1.3932545648232078, 0.053676368673279545], [-2.9970239643019103, -2.7837283183799135, 0.08524408716294365], [-4.47595427893113, -2.6558625932156192, 0.5912091450482926], [-5.036363599701673, -1.5773343835481926, 0.07234477238229159], [-5.17854859712389, -3.730553107169222, 0.2443639867796495], [-4.563626718048291, -2.512914234582642, 2.4603112399854252], [-4.365065377906009, -3.8071850391868187, 3.015169317659496], [-3.808406466198901, -1.3706524210407085, 2.8211950054396913], [-6.087668062292554, -2.1368601745373565, 2.5884672182563873], [-3.015076885801351, -3.287350694771372, -1.1456215214756098], [-2.3510956819310325, -3.610181153269007, 0.8977600360686785], [-2.642782557053619, -0.7109480038655706, 1.1404306775613349], [-2.6677893778920367, -0.7313203152515804, -1.0264488679511405], [-0.3710451618282891, -1.9490799840121056, 1.2553873021032935], [-0.24552532002049696, -2.1655580448419776, -0.8973245689258714], [1.5770424436171642, 0.0, 0.0], [2.292718146893921, 1.3915527243580519, 0.0], [2.341079856722335, 2.0598526928949257, 1.4165023767064706], [1.1453916410070466, 2.0331725558915985, 1.9759511228307933], [3.202217992284651, 1.442225996159817, 2.1984030628363262], [2.721555559149521, 3.3152059779287884, 1.2745358845393997], [1.6292323391939878, 2.2122558673105734, -0.8090479336198866], [3.545586830094389, 1.2600392214310587, -0.42809146886197924], [1.9974224573334811, -0.6906780683055275, 1.053572223549302], [1.9277183224308903, -0.6529932317206316, -1.1102241252095315], [-0.3501493572534735, 0.5705349971623096, -1.1530217920585812], [-0.42668432219275376, 0.7576153073313074, 1.0049834283127264], [-6.6158401798865585, -2.9287936852166743, 2.7612037678207257]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0050', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
