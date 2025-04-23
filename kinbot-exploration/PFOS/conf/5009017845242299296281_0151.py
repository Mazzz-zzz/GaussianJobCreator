import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0151'
logfile = 'conf/5009017845242299296281_0151.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863827, 0.6217394783082105, -1.2501828803165036], [-0.34660204151390067, 2.1278181305643185, -1.566386336981131], [-0.7363023803695165, 3.2026189464980024, -0.4919021834050507], [-2.245266873553545, 3.61521560384561, -0.603146579477167], [-2.640958444907437, 4.178902491939338, 0.5247661900345982], [-2.4192802354782867, 4.467588327805782, -1.6091685125747355], [-3.3507490824222415, 2.1345510510260244, -0.9305681716383417], [-4.68981527198154, 2.5029451046505424, -0.6245966549194926], [-2.9430491354453916, 1.5830044784556052, -2.1694926661774163], [-2.812813595393545, 1.1921710721786434, 0.21096063530076364], [-0.5266683954723089, 2.695152934634475, 0.7195116152968817], [0.004125934750531798, 4.291386672281076, -0.6555498260725737], [0.9809271994195113, 2.1566775839163723, -1.7135621949231834], [-0.9145558447621658, 2.4557522019201006, -2.7205903668258777], [-2.0119541879597165, 0.6042709716797781, -1.0485616399675992], [-0.4099470658637743, -0.1158729665823104, -2.3181096973944295], [1.5770424436171664, 0.0, 0.0], [2.292718146893914, 1.3915527243580557, 0.0], [1.6005215470082423, 2.4407219045638913, -0.9357086002340337], [1.3760692761371616, 1.9145350254105638, -2.125714308102146], [0.4611709133627442, 2.851946778983245, -0.41882897761345406], [2.4048649076934847, 3.4780800111829073, -1.0689303403306902], [3.5400592233304646, 1.2295174124846107, -0.43152105569274335], [2.308846803952242, 1.8960947387583784, 1.2310220414904687], [1.997422457333479, -0.6906780683055241, 1.0535722235493064], [1.9277183224308927, -0.6529932317206322, -1.1102241252095264], [-0.3501493572534774, 0.7132786644586387, 1.0706086973199291], [-0.4266843221927577, -1.249148832966849, 0.15362238828850439], [-3.372354856186968, 1.2632478772648403, 0.9969784798465755]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0151', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
