import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_r13_insertion_ROR_6_5_8_11'
logfile = '5009017845242299296281_r13_insertion_ROR_6_5_8_11.log'

atom = ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H']
geom = [[2.104673, -4.63472, -1.926683], [2.005879, -3.110569, -2.254532], [1.93112, -2.300467, -0.900999], [2.115937, -0.77629, -1.218094], [1.781221, 0.063602, 0.083811], [-0.190454, -0.262778, 0.00653], [1.942962, -0.690643, 1.18406], [1.83046, 1.511908, 0.133779], [2.053798, 2.083298, 1.436489], [2.402385, 2.146202, -1.03673], [-0.129447, 1.725311, -0.202715], [3.356753, -0.530735, -1.61804], [1.35114, -0.389345, -2.229834], [0.771598, -2.539193, -0.281377], [2.8564, -2.714144, -0.04207], [0.94512, -2.841651, -3.012943], [3.043416, -2.71578, -2.988125], [1.903811, -5.457632, -3.252882], [2.274375, -6.957307, -2.984199], [1.833627, -7.826553, -4.209072], [0.520272, -7.916202, -4.317716], [2.27625, -7.356077, -5.357127], [2.282738, -9.063765, -4.116779], [1.7104, -7.406725, -1.869187], [3.581106, -7.090225, -2.765333], [0.65377, -5.340679, -3.689388], [2.645157, -4.962044, -4.241713], [3.275216, -4.917801, -1.35999], [1.204487, -4.994357, -1.015465], [-0.289148, 2.648818, -0.443896]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_r13_insertion_ROR_6_5_8_11', 'label': '5009017845242299296281_r13_insertion_ROR_6_5_8_11', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 18 F\n1 28 F\n1 29 F\n2 3 F\n2 16 F\n2 17 F\n3 4 F\n3 14 F\n3 15 F\n4 5 F\n4 12 F\n4 13 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n8 10 F\n8 11 F\n11 30 F\n18 19 F\n18 26 F\n18 27 F\n19 20 F\n19 24 F\n19 25 F\n20 21 F\n20 22 F\n20 23 F\n'}
Gaussian.command = 'g16 < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy() # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    freq = reader_gauss.read_freq(logfile, ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'])
    zpe = reader_gauss.read_zpe(logfile)
    db.write(mol, name=label, data={'energy': e,'frequencies': np.asarray(freq), 'zpe':zpe, 'status': 'normal'})
except RuntimeError:
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
        db.write(mol, name=label, data={'status': 'error'})

with open(logfile,'a') as f:
    f.write('done\n')
