import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0002'
logfile = 'conf/5009017845242299296281_0002.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863818, 0.6217394783082119, -1.250182880316501], [-0.39761971585595707, -0.07566485901595936, -2.633910119820667], [-0.7491833517666582, -1.5988198371506717, -2.768004365822095], [-0.9314389667275274, -2.0235549081443116, -4.266659819649046], [-0.8151262318136054, -3.335181080650939, -4.378954407094113], [-2.125609869618872, -1.644913545120827, -4.713299062717092], [0.35577392705740385, -1.2298774894092341, -5.377885472781138], [-0.019766177887423762, 0.12416317344822386, -5.5969163398217425], [1.623318271893937, -1.6329111361471018, -4.891827878335688], [0.06381288580124321, -2.0532647478353794, -6.688356686833208], [0.24079307197239347, -2.316267212136122, -2.244356108808552], [-1.879565640212353, -1.8580563759721582, -2.123440869805501], [-1.1305489784219334, 0.5952337900285993, -3.5269594944504687], [0.8878652548597797, 0.09037708117406054, -2.9207161247297284], [-0.256044575953431, 1.8735740976390418, -1.3916701657561963], [-2.007602477187446, 0.6427130616946751, -1.054309216628062], [1.5770424436171646, 0.0, 0.0], [2.2927181468939146, 1.3915527243580554, 0.0], [1.6005215470082494, 2.440721904563896, -0.9357086002340296], [1.3760692761371605, 1.9145350254105762, -2.1257143081021437], [0.46117091336275373, 2.85194677898325, -0.41882897761344795], [2.404864907693494, 3.4780800111829118, -1.0689303403306842], [3.5400592233304637, 1.2295174124846104, -0.4315210556927463], [2.3088468039522487, 1.896094738758375, 1.2310220414904713], [1.997422457333486, -0.6906780683055269, 1.0535722235492981], [1.927718322430894, -0.6529932317206273, -1.1102241252095315], [-0.3501493572534743, 0.7132786644586389, 1.0706086973199314], [-0.4266843221927551, -1.2491488329668516, 0.15362238828850247], [-0.5380780130653403, -1.5633281969790749, -7.26598278448861]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0002', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
