import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0131'
logfile = 'conf/5009017845242299296281_0131.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863838, -1.393559872884596, 0.08664925740765068], [-0.3466020415139057, -2.4204394252486634, -1.0595513872112605], [1.1624535490467036, -2.799902192493109, -1.2603960465750466], [1.322955526797157, -4.17169637078313, -2.0033664407024334], [1.1329590874195934, -5.166849009408311, -1.1547943919415806], [0.44565046019380916, -4.26465604409698, -2.998652510032494], [3.0339402189610625, -4.364524649066058, -2.7496987858344295], [3.2390234210468947, -5.741037082070683, -3.0422486274271114], [3.1906412013888583, -3.32088806457046, -3.693957359137601], [3.8812157953462396, -3.9870416701408957, -1.4768384918147133], [1.7500514572469126, -1.8550381832545337, -1.989060753302593], [1.7633282149005816, -2.8945149144395605, -0.0812319908686853], [-1.0117971004259017, -3.534969835652357, -0.7430628651027298], [-0.8215282967123827, -1.9466956411138756, -2.2051128650057983], [-2.011954187959717, -1.2102165034856942, 0.0009668077395966824], [-0.40994706586377744, -1.949605403411484, 1.2594037813693515], [1.5770424436171688, 0.0, 0.0], [2.29271814689392, 1.391552724358053, 0.0], [2.341079856722326, 2.0598526928949266, 1.416502376706476], [1.1453916410070393, 2.033172555891589, 1.975951122830796], [3.202217992284649, 1.4422259961598232, 2.1984030628363267], [2.721555559149507, 3.315205977928789, 1.274535884539405], [1.6292323391939783, 2.2122558673105743, -0.8090479336198844], [3.545586830094379, 1.2600392214310743, -0.4280914688619865], [1.9974224573334833, -0.6906780683055274, 1.0535722235492964], [1.927718322430895, -0.6529932317206253, -1.1102241252095322], [-0.3501493572534729, 0.570534997162311, -1.1530217920585804], [-0.4266843221927545, 0.7576153073313024, 1.0049834283127308], [3.360896402506907, -3.4319032409464385, -0.8793006766407988]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0131', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
