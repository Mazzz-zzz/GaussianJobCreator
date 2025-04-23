import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0026'
logfile = 'conf/5009017845242299296281_0026.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863816, -1.393559872884597, 0.08664925740765697], [-2.2709622836291916, -1.3932545648232169, 0.05367636867327989], [-2.9970239643019005, -2.783728318379924, 0.08524408716294399], [-4.4759542789311215, -2.6558625932156326, 0.5912091450482931], [-5.036363599701665, -1.5773343835482123, 0.07234477238229249], [-5.178548597123875, -3.7305531071692424, 0.24436398677965082], [-4.5636267180482815, -2.5129142345826576, 2.460311239985427], [-5.8470664672130335, -2.0086793482011642, 2.808105140779786], [-4.008587173203819, -3.7030894512540717, 2.9901893664711485], [-3.5016327596326593, -1.3655548877188255, 2.6507790621022114], [-3.0150768858013435, -3.2873506947713786, -1.145621521475609], [-2.351095681931019, -3.610181153269015, 0.8977600360686805], [-2.6427825570536148, -0.7109480038655811, 1.1404306775613362], [-2.667789377892033, -0.7313203152515892, -1.0264488679511394], [-0.37104516182828207, -1.9490799840121087, 1.255387302103293], [-0.24552532002048844, -2.1655580448419816, -0.897324568925872], [1.577042443617167, 0.0, 0.0], [2.2927181468939177, 1.391552724358057, 0.0], [2.3410798567223274, 2.0598526928949292, 1.4165023767064715], [1.1453916410070408, 2.0331725558915945, 1.975951122830791], [3.2022179922846483, 1.442225996159825, 2.1984030628363236], [2.7215555591495013, 3.3152059779287923, 1.2745358845394064], [1.6292323391939803, 2.212255867310576, -0.8090479336198858], [3.5455868300943862, 1.260039221431073, -0.42809146886197874], [1.997422457333484, -0.6906780683055233, 1.0535722235493017], [1.9277183224308942, -0.6529932317206272, -1.1102241252095326], [-0.35014935725347457, 0.5705349971623072, -1.1530217920585804], [-0.4266843221927563, 0.7576153073313051, 1.004983428312727], [-2.9999597154456636, -1.503982756165492, 3.466330574067868]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0026', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
