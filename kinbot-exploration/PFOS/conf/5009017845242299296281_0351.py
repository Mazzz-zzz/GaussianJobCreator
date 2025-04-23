import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0351'
logfile = 'conf/5009017845242299296281_0351.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863818, 0.6217394783082154, -1.250182880316501], [-2.2709622836291916, 0.6501421835576537, -1.2334320314121774], [-2.9970239643019005, 1.3180406141844456, -2.4534014845326184], [-4.4759542789311215, 0.8159291580463254, -2.5956490472096987], [-5.1688547075294755, 1.6639414173142932, -3.3354234366396467], [-4.503853254513126, -0.38856665038692495, -3.1589531703804945], [-5.326480128927349, 0.6708866980883523, -0.9291461836020357], [-4.864052131593738, -0.5153365265428741, -0.2955328816287247], [-5.286783364380759, 1.955290473341915, -0.3340933450683977], [-6.795257746882602, 0.4289644456798581, -1.4436445437474157], [-3.0150768858013404, 2.6358126881057564, -2.274118452082631], [-2.3510956819310196, 1.0276075788966061, -3.575388609029107], [-2.6427825570536156, -0.6321679360904185, -1.1859143708980966], [-2.667789377892032, 1.2545909529572643, -0.12011753733594861], [-0.37104516182828257, -0.11265730320379994, -2.3156464312139002], [-0.2455253200204879, 1.8598848945507198, -1.426765995739976], [1.5770424436171662, 0.0, 0.0], [2.2927181468939177, 1.3915527243580583, 0.0], [1.6005215470082468, 2.4407219045638984, -0.9357086002340309], [1.376069276137165, 1.9145350254105762, -2.1257143081021432], [0.46117091336274885, 2.8519467789832476, -0.4188289776134523], [2.4048649076934976, 3.47808001118291, -1.0689303403306887], [3.5400592233304655, 1.229517412484611, -0.4315210556927463], [2.308846803952254, 1.896094738758376, 1.2310220414904693], [1.997422457333483, -0.6906780683055269, 1.0535722235493006], [1.9277183224308923, -0.6529932317206324, -1.1102241252095348], [-0.35014935725347346, 0.7132786644586332, 1.0706086973199331], [-0.42668432219275626, -1.2491488329668479, 0.15362238828849978], [-7.432367838459058, 0.8574942754140924, -0.8550739847356815]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0351', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
