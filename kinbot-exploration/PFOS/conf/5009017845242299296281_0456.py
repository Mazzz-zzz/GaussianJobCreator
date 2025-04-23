import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0456'
logfile = 'conf/5009017845242299296281_0456.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863824, 0.771820394576386, 1.1635336229088467], [-0.3976197158559553, 2.3188655045575937, 1.2514273698287424], [-0.7491833517666564, 3.196572017163504, -0.0006164121359367136], [-2.081987459715321, 2.729224887216826, -0.6823825222996559], [-1.8542445281745779, 1.663478425896721, -1.4299248825668842], [-2.9956443258697933, 2.4370379916865277, 0.23889810246250864], [-2.808161535638666, 4.064916675635932, -1.7825016017277127], [-3.4383936490907896, 5.031656144782872, -0.9514460742344191], [-1.8272404471374495, 4.374164343105237, -2.755906178430731], [-3.9213819510009915, 3.200440264992654, -2.4853269513494007], [-0.895313736421607, 4.4602796920374335, 0.3875592899002985], [0.22590414679516663, 3.113269857861801, -0.8965488436803728], [-1.1305489784219331, 2.75681962529853, 2.278967330580892], [0.887865254859781, 2.4842258206717536, 1.5386269105814867], [-0.2560445759534278, 0.2684346684142516, 2.3183978473060174], [-2.0076024771874468, 0.5917020341966386, 1.0837604470856925], [1.5770424436171662, 0.0, 0.0], [2.2927181468939177, 1.3915527243580568, 0.0], [2.341079856722327, 2.059852692894931, 1.416502376706473], [1.1453916410070373, 2.033172555891591, 1.975951122830791], [3.2022179922846474, 1.4422259961598287, 2.1984030628363285], [2.7215555591495084, 3.315205977928793, 1.2745358845394055], [1.629232339193979, 2.212255867310574, -0.8090479336198875], [3.5455868300943836, 1.2600392214310716, -0.4280914688619797], [1.9974224573334813, -0.6906780683055249, 1.0535722235493032], [1.927718322430896, -0.6529932317206286, -1.1102241252095293], [-0.350149357253476, -1.2838136616209446, 0.08241309473865308], [-0.42668432219275676, 0.49153352563554265, -1.1586058166012307], [-4.767867044790215, 3.307472410508246, -2.0292902285792094]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0456', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
