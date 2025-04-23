import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0318'
logfile = 'conf/5009017845242299296281_0318.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863829, -1.3935598728845957, 0.08664925740764695], [-0.3466020415139028, -2.4204394252486625, -1.0595513872112627], [1.1624535490467065, -2.799902192493108, -1.2603960465750477], [1.3229555267971613, -4.171696370783128, -2.0033664407024343], [2.5379668391891173, -4.265983898549034, -2.5145995733771147], [1.1168882370515751, -5.180884981157996, -1.1619541390693056], [0.08826989713528173, -4.354401773287276, -3.40473486582027], [0.026477469084686216, -3.119276365261579, -4.106929945620007], [0.34543438450341707, -5.604174144526686, -4.018976038464472], [-1.2179449472652872, -4.508720715199418, -2.5384029004376494], [1.7500514572469168, -1.8550381832545317, -1.9890607533025926], [1.7633282149005822, -2.8945149144395583, -0.08123199086868653], [-1.0117971004258977, -3.534969835652356, -0.7430628651027362], [-0.8215282967123797, -1.946695641113874, -2.2051128650058027], [-2.011954187959717, -1.2102165034856929, 0.0009668077395913069], [-0.4099470658637792, -1.949605403411485, 1.2594037813693473], [1.5770424436171682, 0.0, 0.0], [2.292718146893921, 1.3915527243580559, 0.0], [3.7823355744197213, 1.3186147352454578, -0.4807937764724422], [4.41959249795878, 0.34925073248438676, 0.14976318527135518], [3.8473329977166433, 1.1104432852891986, -1.779574085222873], [4.368010411051587, 2.468434085335197, -0.20560554420871152], [2.293896732220222, 1.870418904473604, 1.2405689893126353], [1.64946494400089, 2.235214894314339, -0.8029305726284794], [1.997422457333482, -0.6906780683055297, 1.0535722235492988], [1.9277183224308976, -0.6529932317206238, -1.1102241252095317], [-0.3501493572534698, 0.5705349971623122, -1.1530217920585837], [-0.42668432219275404, 0.757615307331305, 1.0049834283127268], [-0.9993205560323984, -4.843566686072109, -1.6574693399210945]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0318', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
