import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0011'
logfile = 'conf/5009017845242299296281_0011.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863824, -1.393559872884601, 0.08664925740765224], [-0.3976197158559557, -2.243200645541638, 1.3824827499919188], [-0.7491833517666568, -1.5977521800128385, 2.7686207779580263], [0.3751315092500395, -0.6189636534602664, 3.2559340987158856], [0.8222268376054804, 0.09417028242405033, 2.237087839494269], [-0.09608992658046996, 0.1968411492370606, 4.194716780097905], [1.8294372022557681, -1.542587331444566, 3.9999620443421753], [2.112852701691653, -2.6651967864903243, 3.1742171237211787], [2.7882050813003083, -0.5689459613893706, 4.371638527694932], [1.128311689997794, -2.048214695844082, 5.316513074207859], [-1.8820373528306351, -0.9113444755835781, 2.648296931279657], [-0.8963077255402097, -2.547907795232933, 3.682944859997029], [-1.1305489784219278, -3.352053415327132, 1.2479921638695695], [0.887865254859781, -2.5746029018458167, 1.3820892141482366], [-0.25604457595342717, -2.142008766053289, -0.9267276815498219], [-2.007602477187446, -1.2344150958913154, -0.029451230457633686], [1.5770424436171644, 0.0, 0.0], [2.292718146893914, 1.3915527243580552, 0.0], [2.3410798567223208, 2.0598526928949292, 1.4165023767064746], [1.1453916410070337, 2.033172555891598, 1.9759511228307916], [3.2022179922846434, 1.4422259961598283, 2.198403062836327], [2.7215555591494973, 3.3152059779288, 1.2745358845394035], [1.629232339193974, 2.212255867310576, -0.8090479336198825], [3.5455868300943774, 1.2600392214310756, -0.42809146886198335], [1.9974224573334831, -0.6906780683055277, 1.053572223549298], [1.9277183224308954, -0.6529932317206268, -1.110224125209533], [-0.3501493572534724, 0.5705349971623093, -1.153021792058582], [-0.4266843221927598, 0.7576153073313014, 1.004983428312726], [0.7947630307666215, -2.9482170764709825, 5.195255523981057]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0011', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
