import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0387'
logfile = 'conf/5009017845242299296281_0387.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863869, -1.393559872884596, 0.08664925740765317], [-0.3466020415139064, -2.4204394252486656, -1.0595513872112574], [-0.7363023803695214, -2.0273092602547957, -2.527598274606108], [0.11232051632789275, -2.815384817700019, -3.5851242808869457], [1.3083470166553106, -2.2653289069216225, -3.701146040093138], [0.2433393664383984, -4.087232389653752, -3.218875929103461], [-0.6978983773361234, -2.8041965493719925, -5.277780033872967], [0.28154973790607757, -3.16907999300268, -6.242130836186471], [-1.9516122589011642, -3.447590983376055, -5.13777333006306], [-0.9663283155263512, -1.2564216221045774, -5.391273525380415], [-2.0209682759934653, -2.3115138840431486, -2.7221562889464326], [-0.5312999584009189, -0.7299509739637089, -2.7153888107570543], [0.980927199419505, -2.562327183726282, -1.0109564779824385], [-0.9145558447621743, -3.583976471922479, -0.7664486088494552], [-2.0119541879597196, -1.2102165034856922, 0.0009668077395966809], [-0.4099470658637807, -1.9496054034114796, 1.2594037813693537], [1.577042443617163, 0.0, 0.0], [2.292718146893916, 1.3915527243580523, 0.0], [1.6005215470082566, 2.4407219045638975, -0.9357086002340294], [1.3760692761371693, 1.9145350254105822, -2.1257143081021406], [0.46117091336276106, 2.851946778983256, -0.41882897761344895], [2.4048649076935025, 3.478080011182909, -1.0689303403306814], [3.5400592233304664, 1.2295174124846018, -0.43152105569274074], [2.3088468039522527, 1.896094738758372, 1.2310220414904716], [1.9974224573334818, -0.6906780683055287, 1.053572223549299], [1.927718322430893, -0.6529932317206266, -1.1102241252095324], [-0.3501493572534765, 0.5705349971623105, -1.1530217920585832], [-0.42668432219276126, 0.7576153073313064, 1.004983428312725], [-1.8000780405866827, -1.092984302833285, -5.853995931718414]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0387', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
