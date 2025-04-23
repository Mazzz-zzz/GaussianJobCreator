import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0165'
logfile = 'conf/5009017845242299296281_0165.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586385, -1.3935598728845984, 0.08664925740765084], [-0.3466020415139042, -2.4204394252486665, -1.059551387211258], [-0.9873929842445245, -3.8498102687678486, -0.9708648395635606], [-1.0274356123761041, -4.556709010365342, -2.3702489598329737], [-1.1749864386042272, -5.860276414715164, -2.2107412199798273], [-2.033642009065594, -4.083012157772827, -3.099673403938122], [0.5505558709515905, -4.266413649690837, -3.343560222087306], [1.6529980492054575, -4.411164957441922, -2.4569378050905235], [0.41632863418113125, -4.971571539097976, -4.564253538020448], [0.3487395885195963, -2.735324150231465, -3.652716756224287], [-0.2614927667443652, -4.595762031236104, -0.14284527771147132], [-2.230877915023897, -3.758011046649702, -0.5172703827719055], [-0.7784534903451411, -1.8458422644301298, -2.1856818928204595], [0.9731582247379432, -2.5497215866483596, -1.121342832061486], [-2.0119541879597183, -1.2102165034856995, 0.0009668077395988365], [-0.4099470658637806, -1.949605403411482, 1.259403781369355], [1.5770424436171635, 0.0, 0.0], [2.292718146893912, 1.3915527243580588, 0.0], [2.3410798567223146, 2.059852692894931, 1.416502376706476], [1.145391641007031, 2.033172555891592, 1.9759511228307918], [3.2022179922846394, 1.4422259961598298, 2.198403062836329], [2.7215555591494915, 3.3152059779287972, 1.2745358845394101], [1.6292323391939765, 2.2122558673105757, -0.8090479336198872], [3.545586830094381, 1.26003922143108, -0.4280914688619783], [1.99742245733348, -0.6906780683055234, 1.0535722235492986], [1.9277183224308951, -0.652993231720628, -1.110224125209531], [-0.35014935725347857, 0.5705349971623082, -1.1530217920585824], [-0.42668432219276314, 0.7576153073313046, 1.0049834283127264], [0.8151274584567616, -2.196106642821906, -2.998737229652932]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0165', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
