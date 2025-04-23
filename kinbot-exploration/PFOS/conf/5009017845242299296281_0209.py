import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0209'
logfile = 'conf/5009017845242299296281_0209.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863844, -1.393559872884597, 0.08664925740765572], [-0.34660204151390467, -2.4204394252486687, -1.0595513872112512], [-0.9873929842445249, -3.849810268767848, -0.9708648395635533], [-1.0274356123761048, -4.556709010365344, -2.370248959832966], [0.09873377709933714, -4.3260521067222895, -3.022201290979427], [-1.1928458587485584, -5.868032928307011, -2.2220292810378575], [-2.4473973794174335, -3.932426367906556, -3.4264890336706015], [-2.5117890942670806, -2.51794247406936, -3.2941257874733245], [-2.3658961109428343, -4.610221076010704, -4.667177813563064], [-3.634643723625211, -4.570482023349106, -2.6117045661381955], [-0.2614927667443673, -4.595762031236105, -0.14284527771146316], [-2.2308779150238984, -3.758011046649701, -0.517270382771897], [-0.7784534903451388, -1.8458422644301333, -2.1856818928204547], [0.9731582247379432, -2.5497215866483636, -1.121342832061477], [-2.011954187959718, -1.210216503485697, 0.0009668077396020592], [-0.4099470658637811, -1.9496054034114783, 1.2594037813693602], [1.5770424436171666, 0.0, 0.0], [2.2927181468939164, 1.3915527243580563, 0.0], [1.6005215470082557, 2.4407219045638953, -0.9357086002340272], [1.3760692761371636, 1.9145350254105764, -2.1257143081021455], [0.4611709133627486, 2.85194677898325, -0.4188289776134485], [2.404864907693497, 3.4780800111829033, -1.0689303403306867], [3.5400592233304655, 1.229517412484609, -0.4315210556927403], [2.3088468039522505, 1.896094738758375, 1.231022041490474], [1.9974224573334793, -0.6906780683055224, 1.0535722235493041], [1.9277183224308916, -0.652993231720624, -1.1102241252095282], [-0.35014935725347607, 0.5705349971623083, -1.1530217920585837], [-0.4266843221927588, 0.7576153073313084, 1.004983428312726], [-3.997120040711748, -3.927237465069321, -1.9865567880162405]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0209', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
