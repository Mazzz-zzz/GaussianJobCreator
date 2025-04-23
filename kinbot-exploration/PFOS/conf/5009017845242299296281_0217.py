import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0217'
logfile = 'conf/5009017845242299296281_0217.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863809, 0.6217394783082131, -1.2501828803165034], [-2.270962283629191, 0.6501421835576504, -1.233432031412179], [-2.9970239643018997, 1.3180406141844418, -2.4534014845326224], [-4.4759542789311215, 0.815929158046319, -2.5956490472097005], [-5.036363599701665, 0.7260147810600387, -1.402184032606565], [-5.178548597123876, 1.653651133263405, -3.352935754365358], [-4.5636267180482815, -0.8742349177524391, -3.4064031846728247], [-3.5410856494537764, -1.6892831007696327, -2.8473528915792836], [-5.9315139537235275, -1.2377206853392295, -3.4524398970509864], [-4.133463285573643, -0.45125614200512937, -4.861254038347501], [-3.015076885801341, 2.635812688105751, -2.2741184520826345], [-2.3510956819310196, 1.0276075788966021, -3.5753886090291083], [-2.6427825570536125, -0.6321679360904238, -1.1859143708981024], [-2.667789377892033, 1.25459095295726, -0.12011753733595047], [-0.3710451618282821, -0.11265730320380458, -2.3156464312139007], [-0.24552532002048846, 1.8598848945507176, -1.4267659957399799], [1.5770424436171666, 0.0, 0.0], [2.2927181468939164, 1.3915527243580537, 0.0], [3.782335574419715, 1.3186147352454562, -0.48079377647244714], [4.419592497958779, 0.3492507324843984, 0.14976318527135168], [3.8473329977166393, 1.1104432852891963, -1.7795740852228783], [4.3680104110515785, 2.468434085335206, -0.20560554420872051], [2.293896732220219, 1.8704189044736057, 1.2405689893126313], [1.6494649440008797, 2.235214894314336, -0.8029305726284859], [1.997422457333484, -0.6906780683055291, 1.0535722235492977], [1.9277183224308936, -0.6529932317206302, -1.1102241252095324], [-0.35014935725347557, 0.7132786644586347, 1.0706086973199331], [-0.426684322192755, -1.2491488329668512, 0.1536223882884979], [-4.6120449169254965, -0.9751396394312181, -5.518877099038965]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0217', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
