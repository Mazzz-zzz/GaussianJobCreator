import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0491'
logfile = 'conf/5009017845242299296281_0491.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863851, 0.7718203945763804, 1.1635336229088473], [-0.3466020415139064, 0.2926212946843425, 2.625937724192379], [1.1624535490467038, 0.3084161010830979, 3.0549844500983028], [1.322955526797157, 0.3508819546540591, 4.614478254324753], [1.1329590874195925, 1.5833432251349522, 5.052019695636848], [0.4456504601938103, -0.4645812287616091, 5.192626727607082], [3.0339402189610625, -0.19904667675479737, 5.154638614451783], [3.0916095975658955, -1.617647645534818, 5.072928505176881], [3.961384938999794, 0.66992129027647, 4.530060568733041], [2.9381651035091867, 0.2160134382784139, 6.670915257409063], [1.7500514572469155, -0.7950580504033825, 2.601040568339857], [1.7633282149005798, 1.3769084895275214, 2.547339442971939], [-1.0117971004259019, 1.1239736000383698, 3.432905111838007], [-0.8215282967123816, -0.9363359387499646, 2.788444311143953], [-2.0119541879597183, 0.6059455318059123, 1.0475948322279915], [-0.4099470658637805, 2.065478369993786, 1.058705916025073], [1.577042443617168, 0.0, 0.0], [2.2927181468939106, 1.3915527243580592, 0.0], [2.34107985672231, 2.059852692894936, 1.416502376706475], [1.1453916410070206, 2.0331725558915856, 1.975951122830794], [3.2022179922846346, 1.4422259961598374, 2.1984030628363316], [2.7215555591494818, 3.3152059779288052, 1.274535884539407], [1.6292323391939671, 2.21225586731058, -0.8090479336198877], [3.5455868300943756, 1.2600392214310925, -0.42809146886197963], [1.9974224573334847, -0.6906780683055233, 1.0535722235493015], [1.927718322430898, -0.6529932317206264, -1.1102241252095282], [-0.35014935725347296, -1.2838136616209481, 0.082413094738651], [-0.4266843221927563, 0.4915335256355456, -1.1586058166012283], [2.673699609851581, -0.5419935707871537, 7.21077242700849]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0491', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
