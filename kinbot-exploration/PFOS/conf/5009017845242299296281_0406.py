import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0406'
logfile = 'conf/5009017845242299296281_0406.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863831, 0.7718203945763845, 1.163533622908849], [-0.3466020415139001, 0.29262129468434045, 2.625937724192381], [-0.7363023803695183, -1.175309686243206, 3.0195004580111457], [-2.2452668735535486, -1.2852675418899022, 3.4324418428267736], [-2.987060025853463, -0.5497865061829501, 2.6229280947079516], [-2.6496830196655092, -2.551038990399625, 3.3761367679707073], [-2.5333983411176346, -0.6864026829068565, 5.187435993583287], [-1.7843602874439564, 0.5072448321049589, 5.378691351324875], [-3.922569383914317, -0.7999498080049874, 5.437771040840226], [-1.8138927449984712, -1.8617787732951645, 5.949896073263335], [-0.5266683954723104, -1.9706918044823205, 1.974315100829184], [0.004125934750528905, -1.5779705333152392, 4.044224788693655], [0.9809271994195118, 0.40564959980989995, 2.724518672905619], [-0.9145558447621631, 1.12822427000237, 3.4870389756753237], [-2.0119541879597143, 0.6059455318059228, 1.0475948322279967], [-0.40994706586377316, 2.0654783699937864, 1.05870591602508], [1.5770424436171657, 0.0, 0.0], [2.2927181468939177, 1.391552724358056, 0.0], [3.7823355744197125, 1.3186147352454591, -0.48079377647244725], [4.419592497958777, 0.34925073248439, 0.14976318527135135], [3.8473329977166326, 1.1104432852892026, -1.7795740852228814], [4.368010411051584, 2.4684340853352036, -0.20560554420871824], [2.2938967322202215, 1.8704189044736017, 1.2405689893126364], [1.6494649440008853, 2.235214894314339, -0.8029305726284831], [1.997422457333482, -0.6906780683055297, 1.0535722235492953], [1.927718322430894, -0.6529932317206208, -1.1102241252095333], [-0.35014935725347907, -1.2838136616209415, 0.08241309473864829], [-0.42668432219275715, 0.49153352563555386, -1.1586058166012265], [-2.27909940989655, -2.062179466497969, 6.774144038255813]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0406', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
