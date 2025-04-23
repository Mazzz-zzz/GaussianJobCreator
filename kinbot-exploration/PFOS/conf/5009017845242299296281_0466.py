import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0466'
logfile = 'conf/5009017845242299296281_0466.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863787, -1.393559872884601, 0.08664925740765224], [-2.27096228362919, -1.39325456482322, 0.05367636867327382], [-3.0203184489305843, -0.657332495932221, 1.2192077454993024], [-2.2823447039971083, 0.6594581441991456, 1.6446984279510897], [-1.2544029753924617, 0.36695518978720537, 2.4220527662875906], [-1.852995246132528, 1.3171334802842678, 0.5713754402676405], [-3.4144078989664113, 1.8184222874767098, 2.5917456117322297], [-4.170179418850305, 1.045383242909505, 3.51567799730869], [-2.6396065708936467, 2.9483758797622, 2.949858098174438], [-4.353747779728866, 2.242435854286523, 1.4007814168970552], [-4.24541884408069, -0.34008637241847073, 0.8099700744271614], [-3.094058093800261, -1.4520145440421888, 2.2791477007209395], [-2.5961406643712612, -0.7866250332999489, -1.0913507571730185], [-2.689777026166647, -2.652128594710344, 0.007528686185010403], [-0.371045161828279, -1.9490799840121122, 1.2553873021032882], [-0.24552532002048677, -2.1655580448419838, -0.8973245689258774], [1.5770424436171673, 0.0, 0.0], [2.292718146893915, 1.3915527243580554, 0.0], [3.7823355744197125, 1.318614735245464, -0.4807937764724409], [4.419592497958776, 0.3492507324843972, 0.14976318527135402], [3.8473329977166375, 1.110443285289203, -1.7795740852228765], [4.368010411051575, 2.4684340853352147, -0.20560554420870997], [2.2938967322202153, 1.8704189044736064, 1.2405689893126333], [1.6494649440008795, 2.23521489431434, -0.8029305726284868], [1.9974224573334842, -0.6906780683055243, 1.0535722235493041], [1.9277183224308967, -0.6529932317206237, -1.110224125209531], [-0.350149357253474, 0.5705349971623063, -1.1530217920585821], [-0.42668432219275826, 0.7576153073313012, 1.0049834283127237], [-4.6177013544088785, 3.1684786301018466, 1.494246247019908]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0466', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
