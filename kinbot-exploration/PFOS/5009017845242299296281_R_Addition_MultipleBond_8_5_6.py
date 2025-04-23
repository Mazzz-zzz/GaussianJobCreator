import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_R_Addition_MultipleBond_8_5_6'
logfile = '5009017845242299296281_R_Addition_MultipleBond_8_5_6.log'

atom = ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H']
geom = [[0.0466700347026793, 0.09964087441494601, 0.022679263176955424], [-0.6911584169457669, -1.2738167495893558, 0.11512222788437114], [-2.2674522170054705, -1.2021724546011872, 0.1176837437485911], [-2.9870695391530666, -0.5855527635887092, -1.136233640806358], [-4.461290502462539, -1.0815107671884927, -1.2515690820817245], [-4.480674075802484, -2.340260101192628, -1.893512523111113], [-5.081857339326028, -1.1480169118157313, -0.09754347459530192], [-5.51295181522398, -0.36557100350629496, -2.6245037833053932], [-5.985969800547627, 0.9218933825446869, -2.287268010507421], [-4.813810837252153, -0.6315150190477452, -3.838657502454507], [-6.738050353036696, -1.3409687432000081, -2.455929385851588], [-3.0323601781247214, 0.7397190708845677, -1.0114498879530902], [-2.329241813969973, -0.9269828880525589, -2.2425704927029773], [-2.6769268721786013, -2.4696064150537045, 0.2112794973353363], [-2.6463653698398844, -0.5187544274156747, 1.1931209917933017], [-0.3684908903993848, -1.8452083006675537, 1.2751634900297313], [-0.3023926821922686, -2.0433633297053455, -0.8983926312804089], [1.5950381061911538, -0.05538979985985355, -0.20204778601593384], [2.380129810483253, 1.251912601368724, 0.16678618740463846], [3.8188745677066036, 1.3067571143371948, -0.44749350363094365], [4.44720475487867, 0.17006740519213442, -0.20658665662607267], [3.7519472493035004, 1.5134058369311052, -1.7465363643423053], [4.478606233594265, 2.299612893240192, 0.11672112454063245], [2.4995284608745973, 1.3121124086310159, 1.488307405370819], [1.7055121438836076, 2.30577726726471, -0.2886904361923839], [2.0526067286830054, -1.0452337226102721, 0.5583571930516941], [1.800407388938162, -0.33233613538329565, -1.4881290757426087], [-0.45163764222205804, 0.7750369647403421, -1.0106705049709825], [-0.16578242739953633, 0.7630705991779829, 1.157285643029482], [-6.771083999215527, -2.1684800346442934, -2.9587405463217697]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'mp2', 'basis': '6-31G', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_R_Addition_MultipleBond_8_5_6', 'label': '5009017845242299296281_R_Addition_MultipleBond_8_5_6', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'guess': 'Mix,Always', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 18 F\n1 28 F\n1 29 F\n2 3 F\n2 16 F\n2 17 F\n3 4 F\n3 14 F\n3 15 F\n4 5 F\n4 12 F\n4 13 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n8 10 F\n8 11 F\n11 30 F\n18 19 F\n18 26 F\n18 27 F\n19 20 F\n19 24 F\n19 25 F\n20 21 F\n20 22 F\n20 23 F\n8 5 6 F\n'}
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
