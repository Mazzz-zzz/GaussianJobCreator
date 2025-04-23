import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0181'
logfile = 'conf/5009017845242299296281_0181.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863842, 0.6217394783082182, -1.2501828803164996], [-0.34660204151390284, 2.1278181305643264, -1.566386336981118], [1.1624535490467065, 2.4914860914100183, -1.7945884035232582], [1.924360652536601, 2.68884256620129, -0.43807332933071386], [1.6581688080331816, 3.8850995556641554, 0.05651039478291353], [1.5600180827406167, 1.7564069008244558, 0.4375511691392213], [3.783231787783966, 2.5518231212816707, -0.6559886921968834], [4.116829961099136, 1.1744706515192669, -0.7744937157620712], [4.161464969515624, 3.543223813537814, -1.5937632382788331], [4.217392380657199, 3.055594975274259, 0.7716881369664406], [1.7462981079693096, 1.5025680769231302, -2.465526734628125], [1.2535978601361863, 3.618508316135836, -2.488809853075984], [-1.0117971004258979, 2.410996235613995, -2.689842246735278], [-0.8215282967123789, 2.8830315798638466, -0.58333144613815], [-2.0119541879597187, 0.6042709716797857, -1.0485616399675868], [-0.40994706586378143, -0.11587296658229987, -2.3181096973944255], [1.577042443617166, 0.0, 0.0], [2.2927181468939173, 1.3915527243580552, 0.0], [1.6005215470082503, 2.4407219045638984, -0.9357086002340262], [1.3760692761371605, 1.9145350254105749, -2.1257143081021415], [0.4611709133627542, 2.851946778983254, -0.4188289776134425], [2.4048649076934954, 3.478080011182911, -1.0689303403306827], [3.5400592233304646, 1.2295174124846038, -0.43152105569274957], [2.308846803952255, 1.8960947387583729, 1.2310220414904702], [1.9974224573334856, -0.6906780683055275, 1.0535722235492984], [1.927718322430891, -0.6529932317206256, -1.1102241252095364], [-0.3501493572534726, 0.7132786644586353, 1.0706086973199362], [-0.4266843221927596, -1.2491488329668483, 0.15362238828850205], [5.034132539369231, 3.5707335766645207, 0.7122733138940291]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0181', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
