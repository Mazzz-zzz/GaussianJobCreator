import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0382'
logfile = 'conf/5009017845242299296281_0382.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863829, 0.7718203945763823, 1.1635336229088502], [-0.34660204151390256, 0.2926212946843428, 2.625937724192381], [-0.9873929842445245, 1.084111519680774, 3.8194659122849357], [-1.027435612376104, 0.22565869267366792, 5.131350240546323], [-2.040447113782524, -0.6216591292156103, 5.082748287868308], [0.10871386352754484, -0.45017867605863776, 5.2773237312282335], [-1.2388794809575765, 1.2954728144902563, 6.658587146306556], [-1.6618516099518503, 0.46196952615692655, 7.730407025876611], [-0.11466179393237848, 2.1545409167303244, 6.717691795373941], [-2.467331563379808, 2.1508259644365064, 6.168873862967697], [-0.2614927667443652, 2.1741733763092714, 4.051469307654177], [-2.230877915023899, 1.431036231219077, 3.513168225487142], [-0.778453490345139, -0.9699349115591109, 2.6913872387857145], [0.9731582247379446, 0.30374941440734815, 2.7687950826457897], [-2.0119541879597183, 0.6059455318059167, 1.0475948322279949], [-0.40994706586377944, 2.065478369993784, 1.0587059160250756], [1.577042443617165, 0.0, 0.0], [2.292718146893913, 1.3915527243580583, 0.0], [3.7823355744197156, 1.3186147352454665, -0.48079377647244276], [4.419592497958774, 0.34925073248439564, 0.14976318527135185], [3.8473329977166353, 1.110443285289197, -1.7795740852228772], [4.368010411051577, 2.4684340853352094, -0.20560554420871974], [2.2938967322202144, 1.8704189044736075, 1.2405689893126348], [1.6494649440008788, 2.235214894314339, -0.8029305726284873], [1.9974224573334833, -0.690678068305524, 1.0535722235493004], [1.9277183224308987, -0.6529932317206227, -1.1102241252095286], [-0.3501493572534761, -1.2838136616209437, 0.08241309473865072], [-0.42668432219275904, 0.49153352563554653, -1.1586058166012254], [-2.495056719118269, 2.1763241099278705, 5.202155937711435]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0382', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
