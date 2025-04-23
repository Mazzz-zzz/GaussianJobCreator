import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0293'
logfile = 'conf/5009017845242299296281_0293.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863851, 0.6217394783082199, -1.250182880316496], [-0.34660204151390633, 2.12781813056433, -1.5663863369811133], [1.1624535490467038, 2.491486091410022, -1.794588403523254], [1.322955526797157, 3.8208144161290867, -2.611111813622317], [1.1329590874195923, 3.583505784273377, -3.897225303695263], [0.44565046019381027, 4.7292372728586045, -2.1939742175745796], [3.033940218961062, 4.563571325820874, -2.4049398286173482], [3.0916095975658955, 5.202108779832811, -1.1355402971831812], [3.961384938999794, 3.5881868880667724, -2.845199140281977], [2.9381651035091862, 5.6691753602702635, -3.522530753812439], [1.7500514572469152, 2.6500962336579246, -0.6119798150372594], [1.7633282149005796, 1.51760642491205, -2.466107452103256], [-1.0117971004259019, 2.4109962356139976, -2.6898422467352736], [-0.8215282967123816, 2.8830315798638475, -0.5833314461381448], [-2.0119541879597183, 0.6042709716797858, -1.0485616399675848], [-0.4099470658637805, -0.11587296658229543, -2.3181096973944295], [1.577042443617168, 0.0, 0.0], [2.292718146893916, 1.3915527243580563, 0.0], [2.341079856722326, 2.0598526928949275, 1.4165023767064742], [1.1453916410070395, 2.033172555891586, 1.9759511228307967], [3.2022179922846528, 1.4422259961598294, 2.1984030628363227], [2.7215555591494978, 3.3152059779287963, 1.274535884539405], [1.6292323391939678, 2.2122558673105743, -0.8090479336198838], [3.545586830094379, 1.2600392214310803, -0.4280914688619871], [1.997422457333486, -0.6906780683055286, 1.053572223549295], [1.927718322430892, -0.6529932317206201, -1.1102241252095384], [-0.35014935725347374, 0.7132786644586343, 1.070608697319939], [-0.4266843221927563, -1.2491488329668514, 0.15362238828849792], [2.673699609851581, 6.515708888091315, -3.1360060125147]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0293', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
