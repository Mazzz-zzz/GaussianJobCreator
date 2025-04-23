import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0070'
logfile = 'conf/5009017845242299296281_0070.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863811, -1.3935598728845997, 0.08664925740765092], [-0.3466020415139006, -2.4204394252486643, -1.0595513872112623], [1.1624535490467085, -2.799902192493107, -1.2603960465750472], [1.9172819807099781, -2.9259596001371135, 0.10851911797418273], [3.0265100569963668, -3.6259410075683345, -0.05316281393412801], [2.2192278579630744, -1.7209079512393475, 0.5834709589556054], [0.8705838155292752, -3.7890645027347674, 1.4051054858009635], [1.728184502501133, -4.219182035817121, 2.454844013840466], [-0.2831861675103228, -2.99139871374414, 1.5997869871275574], [0.4444328891869027, -5.044365938219233, 0.554719500581781], [1.2315520144708556, -3.970118415502822, -1.8887358797759235], [1.768541893208143, -1.8669832707228657, -1.983409881825591], [-1.0117971004258943, -3.5349698356523604, -0.7430628651027306], [-0.8215282967123791, -1.9466956411138783, -2.2051128650057996], [-2.011954187959716, -1.2102165034856984, 0.0009668077395966857], [-0.4099470658637744, -1.9496054034114847, 1.259403781369352], [1.577042443617166, 0.0, 0.0], [2.292718146893918, 1.3915527243580557, 0.0], [2.34107985672233, 2.0598526928949297, 1.4165023767064704], [1.145391641007047, 2.033172555891596, 1.9759511228307947], [3.2022179922846554, 1.4422259961598265, 2.198403062836325], [2.7215555591495084, 3.315205977928793, 1.274535884539402], [1.6292323391939783, 2.2122558673105743, -0.8090479336198907], [3.545586830094383, 1.2600392214310734, -0.4280914688619894], [1.9974224573334847, -0.6906780683055261, 1.053572223549295], [1.927718322430895, -0.6529932317206221, -1.1102241252095348], [-0.35014935725347474, 0.57053499716231, -1.1530217920585797], [-0.42668432219275326, 0.7576153073313026, 1.0049834283127292], [-0.4596840006967661, -5.308597096397844, 0.7754386233841987]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0070', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
