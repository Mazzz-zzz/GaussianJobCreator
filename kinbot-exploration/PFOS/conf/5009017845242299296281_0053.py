import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0053'
logfile = 'conf/5009017845242299296281_0053.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863809, -1.393559872884599, 0.08664925740765461], [-0.346602041513899, -2.4204394252486687, -1.0595513872112565], [-0.7363023803695162, -2.0273092602548006, -2.527598274606105], [-2.245266873553545, -2.3299480619557147, -2.8292952633496196], [-2.9870600258534603, -1.9966291092255306, -1.7875931283662994], [-2.6496830196655052, -1.6483007125135147, -3.897332955716036], [-2.5333983411176324, -4.149250009495477, -3.1881601574147673], [-3.9242325103754014, -4.41411814827757, -3.054533099925813], [-1.7723615231531522, -4.466813142699635, -4.339391406895107], [-1.782885961350609, -4.748951130016888, -1.94010873072162], [-0.5266683954723086, -0.7244611301521547, -2.6938267161260767], [0.004125934750531798, -2.7134161389658473, -3.388674962621098], [0.980927199419513, -2.5623271837262833, -1.010956477982439], [-0.9145558447621623, -3.583976471922483, -0.7664486088494511], [-2.0119541879597147, -1.2102165034857035, 0.0009668077396031392], [-0.40994706586377255, -1.9496054034114843, 1.2594037813693566], [1.5770424436171673, 0.0, 0.0], [2.2927181468939093, 1.39155272435806, 0.0], [2.341079856722307, 2.0598526928949292, 1.416502376706483], [1.145391641007024, 2.0331725558915843, 1.9759511228307984], [3.2022179922846354, 1.4422259961598325, 2.198403062836336], [2.72155555914948, 3.315205977928801, 1.274535884539416], [1.6292323391939616, 2.2122558673105805, -0.8090479336198829], [3.5455868300943747, 1.2600392214310943, -0.42809146886197685], [1.9974224573334856, -0.6906780683055224, 1.053572223549299], [1.9277183224308958, -0.6529932317206242, -1.110224125209531], [-0.35014935725347873, 0.5705349971623066, -1.1530217920585777], [-0.426684322192758, 0.7576153073313019, 1.0049834283127301], [-2.410083600011605, -4.931097619741068, -1.2263803755084985]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0053', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
