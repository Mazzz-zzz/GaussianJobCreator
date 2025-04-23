import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0192'
logfile = 'conf/5009017845242299296281_0192.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863829, 0.6217394783082193, -1.250182880316496], [-0.3466020415139028, 2.1278181305643273, -1.566386336981116], [1.1624535490467065, 2.49148609141002, -1.7945884035232547], [1.3229555267971613, 3.820814416129083, -2.6111118136223195], [2.537966839189117, 4.31069906016461, -2.437150641590282], [1.1168882370515751, 3.5967242930454937, -3.9058009382334338], [0.08826989713528173, 5.125787773594594, -2.068655121040655], [0.026477469084686216, 5.116343847100759, -0.6479076009309124], [0.34543438450341707, 6.282622418774521, -2.8438691571597996], [-1.2179449472652872, 4.452681754418814, -2.6354652277130155], [1.7500514572469168, 2.6500962336579232, -0.6119798150372616], [1.7633282149005822, 1.517606424912047, -2.4661074521032558], [-1.0117971004258977, 2.4109962356139967, -2.689842246735275], [-0.8215282967123797, 2.8830315798638475, -0.5833314461381488], [-2.011954187959717, 0.6042709716797849, -1.0485616399675843], [-0.4099470658637792, -0.1158729665822964, -2.318109697394428], [1.5770424436171682, 0.0, 0.0], [2.2927181468939195, 1.3915527243580552, 0.0], [1.600521547008255, 2.4407219045639006, -0.9357086002340254], [1.376069276137165, 1.914535025410577, -2.125714308102141], [0.46117091336276195, 2.851946778983254, -0.41882897761344196], [2.4048649076934994, 3.4780800111829113, -1.0689303403306836], [3.5400592233304673, 1.229517412484603, -0.4315210556927489], [2.3088468039522585, 1.896094738758371, 1.2310220414904702], [1.997422457333488, -0.6906780683055295, 1.0535722235492966], [1.9277183224308938, -0.6529932317206237, -1.1102241252095357], [-0.3501493572534703, 0.7132786644586351, 1.070608697319938], [-0.42668432219275404, -1.2491488329668512, 0.15362238828850242], [-0.9993205560323984, 3.8571938974015483, -3.365917125101906]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0192', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
