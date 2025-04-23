import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0189'
logfile = 'conf/5009017845242299296281_0189.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863828, -1.3935598728845968, 0.08664925740765198], [-0.3976197158559559, -2.243200645541637, 1.3824827499919152], [1.0879414097563267, -2.676195207864829, 1.6417185410974604], [1.473042089945535, -3.9541613281908345, 0.8182592400881694], [0.9614495666061577, -3.875818325301181, -0.39773311475559364], [2.7953494776680783, -4.070811596542063, 0.7358218224855981], [0.8248630198382829, -5.522044717330672, 1.6202997018164225], [1.6664714540296413, -5.838737038135043, 2.722040245343579], [-0.5819280292581736, -5.3868141226324875, 1.7097893776883573], [1.1337571966355793, -6.523456817456966, 0.44458680188278227], [1.894753649284141, -1.6828188678150051, 1.279435489212164], [1.2598103705638592, -2.9456573866898466, 2.929461138040974], [-0.7843720924895148, -1.4713296274432146, 2.4020244202192775], [-1.158198376361053, -3.330856542149739, 1.3537778366861803], [-0.25604457595342933, -2.142008766053283, -0.9267276815498251], [-2.007602477187446, -1.2344150958913067, -0.02945123045763403], [1.5770424436171653, 0.0, 0.0], [2.2927181468939173, 1.3915527243580559, 0.0], [2.3410798567223217, 2.0598526928949332, 1.4165023767064737], [1.1453916410070397, 2.0331725558915967, 1.975951122830796], [3.2022179922846497, 1.4422259961598294, 2.1984030628363236], [2.721555559149501, 3.3152059779287986, 1.2745358845394046], [1.6292323391939738, 2.212255867310578, -0.8090479336198799], [3.5455868300943783, 1.260039221431081, -0.4280914688619834], [1.9974224573334856, -0.6906780683055271, 1.053572223549297], [1.927718322430892, -0.6529932317206251, -1.1102241252095366], [-0.35014935725347807, 0.5705349971623109, -1.153021792058579], [-0.4266843221927557, 0.7576153073313041, 1.0049834283127295], [1.9831461677131388, -6.962635181557486, 0.5916348450992769]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0189', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
