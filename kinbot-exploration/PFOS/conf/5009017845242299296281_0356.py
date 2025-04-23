import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0356'
logfile = 'conf/5009017845242299296281_0356.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863829, -1.393559872884596, 0.0866492574076457], [-0.3976197158559581, -2.2432006455416387, 1.382482749991911], [1.0879414097563231, -2.676195207864833, 1.6417185410974564], [1.4730420899455323, -3.954161328190836, 0.818259240088166], [1.0004737759019127, -5.033376622340047, 1.416987459327409], [0.9773040172258539, -3.878060131557083, -0.41354935367959517], [3.3329235748106685, -4.146760933982903, 0.6588553969111792], [3.6096714118651945, -5.486548726400827, 0.2705123606873445], [3.8065427668236027, -3.006068515267075, -0.03386564253344726], [3.7106516445294835, -3.972982230044219, 2.1779554914450765], [1.894753649284139, -1.6828188678150073, 1.2794354892121609], [1.2598103705638528, -2.9456573866898506, 2.9294611380409674], [-0.7843720924895168, -1.471329627443218, 2.4020244202192735], [-1.1581983763610575, -3.3308565421497387, 1.3537778366861735], [-0.25604457595342994, -2.1420087660532823, -0.9267276815498292], [-2.0076024771874463, -1.2344150958913058, -0.02945123045763895], [1.5770424436171657, 0.0, 0.0], [2.2927181468939186, 1.3915527243580539, 0.0], [2.3410798567223248, 2.0598526928949314, 1.4165023767064708], [1.1453916410070428, 2.033172555891597, 1.9759511228307902], [3.20221799228465, 1.4422259961598316, 2.198403062836322], [2.721555559149504, 3.3152059779287955, 1.2745358845394001], [1.6292323391939725, 2.212255867310575, -0.8090479336198849], [3.5455868300943765, 1.2600392214310776, -0.42809146886198834], [1.997422457333484, -0.6906780683055294, 1.053572223549296], [1.9277183224308954, -0.6529932317206238, -1.1102241252095362], [-0.35014935725347274, 0.5705349971623127, -1.15302179205858], [-0.426684322192757, 0.7576153073313023, 1.004983428312727], [4.55195139009703, -3.503011124593535, 2.2634206003810515]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0356', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
