import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0182'
logfile = 'conf/5009017845242299296281_0182.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863858, -1.3935598728845966, 0.08664925740765321], [-2.2709622836291956, -1.393254564823212, 0.05367636867327847], [-3.020318448930586, -0.657332495932207, 1.2192077454993075], [-3.131423238500186, -1.5575296086444679, 2.4986162193657293], [-1.9982160731709866, -2.2106058741423933, 2.687917183353885], [-3.396984742281267, -0.8121282410883275, 3.5676129123016715], [-4.510308809776061, -2.8211047900622277, 2.344617279930162], [-4.300219575237264, -3.830322790937119, 3.3242970990807823], [-5.717986111385332, -2.0961031467770477, 2.198860893930951], [-4.129514471083674, -3.402243715468467, 0.9311761418538339], [-2.3454389866601115, 0.4439741098051406, 1.5364663530591338], [-4.247573342377323, -0.33444956237991935, 0.8318812175920345], [-2.596140664371266, -0.7866250332999382, -1.0913507571730148], [-2.6897770261666554, -2.6521285947103332, 0.0075286861850174394], [-0.37104516182828695, -1.9490799840121065, 1.2553873021032915], [-0.24552532002049418, -2.165558044841981, -0.8973245689258718], [1.5770424436171648, 0.0, 0.0], [2.2927181468939177, 1.3915527243580552, 0.0], [3.782335574419718, 1.3186147352454534, -0.4807937764724495], [4.419592497958775, 0.34925073248438565, 0.1497631852713427], [3.8473329977166384, 1.1104432852891892, -1.7795740852228823], [4.368010411051586, 2.4684340853351996, -0.20560554420872168], [2.2938967322202197, 1.870418904473606, 1.2405689893126273], [1.6494649440008857, 2.2352148943143355, -0.8029305726284902], [1.99742245733348, -0.6906780683055233, 1.0535722235493021], [1.9277183224308916, -0.6529932317206321, -1.1102241252095348], [-0.35014935725347623, 0.5705349971623079, -1.1530217920585843], [-0.42668432219275654, 0.7576153073313094, 1.0049834283127217], [-3.593604541986497, -4.201645320883573, 1.0297662420201967]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0182', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
