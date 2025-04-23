import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_intra_H_migration_9_30'
logfile = '5009017845242299296281_intra_H_migration_9_30.log'

atom = ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H']
geom = [[2.453041, -4.015472, -5.330742], [2.473899, -3.728366, -3.792508], [1.930039, -2.301938, -3.541996], [1.980734, -1.96136, -1.961314], [1.432224, -0.597157, -1.730947], [2.103162, 0.374854, -2.448219], [0.167514, -0.453205, -2.0987], [1.58668, 0.136891, 0.302082], [0.704042, -0.326648, 0.879715], [3.033935, -0.108361, 0.624891], [1.104034, 1.337809, 0.207874], [1.316334, -2.887734, -1.291162], [3.22614, -2.096736, -1.522472], [2.657233, -1.375269, -4.229817], [0.686313, -2.149499, -4.008809], [3.723258, -3.849067, -3.328052], [1.798017, -4.627626, -3.124547], [2.757556, -5.552109, -5.570161], [2.892122, -5.800368, -7.149677], [3.357258, -7.256698, -7.389487], [4.602494, -7.451699, -7.033517], [2.641973, -8.137655, -6.731351], [3.279714, -7.595498, -8.649002], [3.72869, -4.938187, -7.695329], [1.733447, -5.563067, -7.748559], [3.852777, -5.931423, -4.943774], [1.796067, -6.308122, -5.073337], [1.285929, -3.69199, -5.877445], [3.3391, -3.255973, -5.962192], [0.224869, 0.773728, 0.802945]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_intra_H_migration_9_30', 'label': '5009017845242299296281_intra_H_migration_9_30', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': ''}
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
