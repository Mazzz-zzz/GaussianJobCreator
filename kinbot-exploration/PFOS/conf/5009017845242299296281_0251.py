import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0251'
logfile = 'conf/5009017845242299296281_0251.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863844, -1.3935598728845984, 0.08664925740765457], [-0.3466020415139049, -2.42043942524867, -1.0595513872112508], [-0.9873929842445254, -3.84981026876785, -0.9708648395635536], [-1.0274356123761064, -4.5567090103653465, -2.3702489598329652], [-2.0404471137825277, -4.090959573728019, -3.079746742329375], [0.10871386352754246, -4.345207077208822, -3.028528035322929], [-1.238879480957578, -6.414242029259143, -2.207381205892571], [-2.2191735832519255, -6.6643637433100915, -1.2079345545013445], [-1.2941925262646914, -6.938572902463455, -3.521677357477773], [0.18152008934522704, -6.754064281666181, -1.6178657037540534], [-0.2614927667443673, -4.595762031236105, -0.14284527771146316], [-2.230877915023899, -3.758011046649702, -0.5172703827718971], [-0.7784534903451408, -1.8458422644301349, -2.185681892820454], [0.9731582247379439, -2.549721586648366, -1.1213428320614778], [-2.011954187959721, -1.2102165034856966, 0.0009668077396031337], [-0.40994706586378105, -1.9496054034114778, 1.2594037813693597], [1.5770424436171633, 0.0, 0.0], [2.2927181468939177, 1.3915527243580539, 0.0], [1.600521547008256, 2.4407219045639, -0.9357086002340305], [1.3760692761371716, 1.9145350254105828, -2.1257143081021415], [0.4611709133627573, 2.8519467789832538, -0.41882897761344806], [2.404864907693505, 3.478080011182908, -1.0689303403306785], [3.5400592233304673, 1.2295174124846024, -0.4315210556927357], [2.3088468039522487, 1.8960947387583702, 1.231022041490471], [1.9974224573334811, -0.6906780683055241, 1.0535722235493004], [1.9277183224308942, -0.6529932317206291, -1.110224125209528], [-0.3501493572534766, 0.5705349971623055, -1.1530217920585832], [-0.4266843221927601, 0.7576153073313069, 1.0049834283127237], [0.8206069554527772, -6.071630932364978, -1.8664921062666733]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0251', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
