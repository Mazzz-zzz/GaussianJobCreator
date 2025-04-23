import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0145'
logfile = 'conf/5009017845242299296281_0145.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586381, -1.393559872884598, 0.0866492574076533], [-0.3466020415139008, -2.420439425248669, -1.059551387211254], [-0.9873929842445209, -3.8498102687678504, -0.9708648395635574], [-0.18477695124074436, -4.7876250073543885, -0.0034475673086363707], [0.8897635883676692, -5.251777005738081, -0.6169844400485593], [0.18128897051813653, -4.122301251114864, 1.0884499008005055], [-1.2119219398679437, -6.254814906871081, 0.5569008084445375], [-1.9175735786971912, -6.759974027311122, -0.5697318004674498], [-0.382611162427967, -7.042364963812732, 1.3919848245058652], [-2.225283750583378, -5.482697679467235, 1.4829005155388462], [-2.230506898554333, -3.7359427839966646, -0.5118511747699833], [-1.003631458776226, -4.4109250696587035, -2.1730779211988818], [-0.7784534903451346, -1.845842264430132, -2.185681892820458], [0.9731582247379498, -2.5497215866483645, -1.1213428320614773], [-2.011954187959715, -1.2102165034856978, 0.0009668077395988351], [-0.40994706586377727, -1.9496054034114811, 1.2594037813693546], [1.5770424436171675, 0.0, 0.0], [2.292718146893918, 1.391552724358054, 0.0], [1.6005215470082517, 2.4407219045638984, -0.9357086002340265], [1.3760692761371691, 1.9145350254105755, -2.125714308102145], [0.4611709133627464, 2.851946778983254, -0.4188289776134436], [2.4048649076935, 3.478080011182904, -1.0689303403306831], [3.540059223330468, 1.2295174124846047, -0.43152105569273935], [2.3088468039522514, 1.8960947387583742, 1.2310220414904691], [1.9974224573334836, -0.6906780683055251, 1.0535722235493017], [1.9277183224308976, -0.6529932317206286, -1.1102241252095308], [-0.35014935725347207, 0.5705349971623055, -1.1530217920585832], [-0.42668432219275443, 0.7576153073313069, 1.0049834283127237], [-3.035262155625421, -5.276125657886656, 0.9958426442733639]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0145', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
