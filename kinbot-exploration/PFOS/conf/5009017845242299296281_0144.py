import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0144'
logfile = 'conf/5009017845242299296281_0144.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863844, 0.7718203945763833, 1.1635336229088493], [-0.39761971585595723, 2.318865504557592, 1.251427369828748], [-0.7491833517666582, 3.196572017163505, -0.0006164121359281963], [0.3751315092500395, 3.1292034692660797, -1.0919288014421265], [-0.10815375808507399, 3.517453645863822, -2.25906165327419], [1.3965384052652348, 3.910843918740281, -0.7533048196106786], [1.0485429431935351, 1.3891449930549475, -1.292868770375519], [1.9279419990532713, 1.1244117256552761, -0.20704898537669397], [-0.04579798938659392, 0.5630900429513476, -1.6469985690891982], [1.907047804188475, 1.6260551230508153, -2.5918628518718134], [-1.8820373528306367, 2.7491646570443424, -0.5349009981858484], [-0.8963077255402104, 4.463477707111217, 0.36508044717360233], [-1.1305489784219285, 2.756819625298525, 2.2789673305809], [0.8878652548597809, 2.484225820671752, 1.538626910581493], [-0.2560445759534276, 0.26843466841424324, 2.318397847306017], [-2.007602477187446, 0.5917020341966366, 1.0837604470856959], [1.577042443617165, 0.0, 0.0], [2.292718146893915, 1.3915527243580539, 0.0], [3.782335574419713, 1.3186147352454638, -0.48079377647244437], [4.419592497958775, 0.34925073248439376, 0.14976318527135213], [3.8473329977166344, 1.1104432852891983, -1.7795740852228779], [4.3680104110515785, 2.468434085335207, -0.20560554420871635], [2.293896732220216, 1.8704189044736057, 1.2405689893126344], [1.6494649440008788, 2.2352148943143395, -0.802930572628487], [1.9974224573334816, -0.690678068305528, 1.0535722235492995], [1.927718322430894, -0.6529932317206244, -1.110224125209532], [-0.3501493572534785, -1.283813661620946, 0.08241309473865085], [-0.42668432219276053, 0.4915335256355454, -1.1586058166012267], [1.882085252168938, 0.842470812958887, -3.158726933173875]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0144', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
