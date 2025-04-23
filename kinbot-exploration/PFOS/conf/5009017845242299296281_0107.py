import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0107'
logfile = 'conf/5009017845242299296281_0107.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863838, -1.393559872884596, 0.0866492574076482], [-2.2709622836291934, -1.3932545648232144, 0.053676368673272364], [-3.020318448930585, -0.6573324959322158, 1.2192077454993058], [-3.131423238500183, -1.5575296086444825, 2.498616219365724], [-3.378446400292141, -0.8053658014461169, 3.556789696521602], [-4.105859700541783, -2.450930597679871, 2.3529855399487385], [-1.5468226482012546, -2.5097665945431706, 2.8210274831668882], [-1.4940447140016708, -3.605063865628989, 1.9155058735828951], [-0.5142916547324224, -1.5519140929025792, 2.9685235993026833], [-1.8980797259386455, -3.0569121857670987, 4.255528370985577], [-2.345438986660113, 0.443974109805129, 1.5364663530591347], [-4.247573342377322, -0.3344495623799274, 0.8318812175920354], [-2.5961406643712666, -0.7866250332999394, -1.0913507571730194], [-2.6897770261666514, -2.652128594710337, 0.007528686185005672], [-0.3710451618282838, -1.9490799840121131, 1.2553873021032815], [-0.24552532002049143, -2.1655580448419753, -0.897324568925882], [1.577042443617165, 0.0, 0.0], [2.292718146893917, 1.3915527243580543, 0.0], [1.6005215470082585, 2.440721904563901, -0.9357086002340282], [1.3760692761371693, 1.9145350254105806, -2.1257143081021415], [0.46117091336276217, 2.8519467789832573, -0.4188289776134497], [2.404864907693509, 3.478080011182908, -1.068930340330687], [3.5400592233304615, 1.2295174124846016, -0.43152105569275023], [2.3088468039522554, 1.8960947387583724, 1.2310220414904716], [1.997422457333482, -0.69067806830553, 1.0535722235492968], [1.9277183224308931, -0.6529932317206211, -1.1102241252095353], [-0.3501493572534762, 0.570534997162314, -1.153021792058581], [-0.42668432219275615, 0.7576153073313036, 1.0049834283127268], [-2.5526952281901254, -2.4862346376331423, 4.681871034771242]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0107', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
