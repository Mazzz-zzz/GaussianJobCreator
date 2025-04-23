import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0428'
logfile = 'conf/5009017845242299296281_0428.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863863, -1.3935598728845997, 0.08664925740764844], [-0.3976197158559574, -2.2432006455416413, 1.3824827499919126], [1.0879414097563262, -2.676195207864833, 1.6417185410974564], [1.473042089945535, -3.954161328190836, 0.8182592400881623], [2.7879690257787106, -4.0509445055325335, 0.7282475015146945], [0.9875700214111528, -5.045551894983047, 1.403138323269409], [0.7750378848301898, -3.898398531812576, -0.9228162813342629], [-0.6091697343464534, -4.217985433125179, -0.8575132390411538], [1.2916422504512648, -2.7309923054974377, -1.5357105657177088], [1.5508075901178944, -5.132826845412615, -1.5185601296024578], [1.8947536492841406, -1.6828188678150071, 1.2794354892121607], [1.2598103705638566, -2.9456573866898528, 2.9294611380409648], [-0.7843720924895178, -1.4713296274432217, 2.402024420219275], [-1.1581983763610542, -3.330856542149742, 1.3537778366861748], [-0.2560445759534314, -2.142008766053281, -0.9267276815498285], [-2.0076024771874486, -1.234415095891311, -0.02945123045763413], [1.5770424436171633, 0.0, 0.0], [2.2927181468939084, 1.3915527243580563, 0.0], [3.7823355744197125, 1.318614735245463, -0.48079377647244104], [4.419592497958774, 0.3492507324843942, 0.14976318527134935], [3.8473329977166344, 1.1104432852892054, -1.779574085222874], [4.36801041105158, 2.468434085335204, -0.2056055442087119], [2.293896732220221, 1.8704189044735997, 1.2405689893126377], [1.6494649440008806, 2.235214894314341, -0.8029305726284761], [1.9974224573334818, -0.6906780683055285, 1.0535722235492964], [1.9277183224308918, -0.6529932317206236, -1.110224125209538], [-0.35014935725347546, 0.5705349971623134, -1.1530217920585784], [-0.4266843221927618, 0.7576153073312999, 1.0049834283127295], [2.3658318764749176, -5.289345721564876, -1.0213683750929288]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0428', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
