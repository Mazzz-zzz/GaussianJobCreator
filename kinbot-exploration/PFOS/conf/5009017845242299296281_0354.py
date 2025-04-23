import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0354'
logfile = 'conf/5009017845242299296281_0354.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863832, 0.6217394783082135, -1.2501828803165], [-0.39761971585595646, -0.07566485901595409, -2.633910119820667], [-1.1233882121466352, 0.4940238789133692, -3.9028837901505367], [-2.5854106473305376, -0.058297607237880016, -4.033331227674649], [-2.560441758543041, -1.2811201495928315, -4.53391441514531], [-3.1805547950289785, -0.07857504986565969, -2.8440910041143335], [-3.636080398918967, 1.009706014139254, -5.163430634336482], [-4.018635167410041, 2.175579655243872, -4.444527643737474], [-2.9923978927211357, 1.0309022274710917, -6.42461274783289], [-4.884444232150147, 0.05668244722504319, -5.281228981504345], [-1.1808982350936241, 1.819319269574914, -3.8065703577649517], [-0.45262953286663543, 0.15318954183902528, -4.995744508783357], [0.9186702629217401, 0.05926285465343761, -2.8183083047232302], [-0.6863393184328429, -1.3666682794365124, -2.5232286857846646], [-0.2560445759534289, 1.873574097639044, -1.3916701657561925], [-2.0076024771874463, 0.6427130616946772, -1.0543092166280603], [1.5770424436171657, 0.0, 0.0], [2.292718146893915, 1.3915527243580557, 0.0], [1.6005215470082472, 2.440721904563902, -0.9357086002340265], [1.3760692761371618, 1.9145350254105782, -2.1257143081021432], [0.4611709133627506, 2.851946778983251, -0.41882897761344995], [2.404864907693493, 3.4780800111829153, -1.0689303403306842], [3.540059223330467, 1.2295174124846135, -0.43152105569274674], [2.308846803952252, 1.8960947387583753, 1.2310220414904702], [1.9974224573334864, -0.6906780683055279, 1.0535722235492977], [1.9277183224308954, -0.652993231720622, -1.1102241252095344], [-0.3501493572534721, 0.713278664458636, 1.0706086973199331], [-0.42668432219275654, -1.2491488329668472, 0.15362238828850416], [-5.562097152855033, 0.3095495684656284, -4.638728300764137]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0354', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
