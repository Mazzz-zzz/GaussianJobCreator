import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_R_Addition_MultipleBond_8_5_4'
logfile = '5009017845242299296281_R_Addition_MultipleBond_8_5_4.log'

atom = ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H']
geom = [[-0.42399676864810254, 1.3965641753950935, -1.2839757557178635], [-0.6659113194279795, -0.06951831556247076, -0.7960609845666181], [-2.133774520025864, -0.4605606675761215, -0.3985958618651956], [-3.2797083789237216, -0.3961725579775133, -1.5058539890507885], [-4.3722889856365805, -1.2857057150684383, -1.2507606850714739], [-4.264965054817134, -2.569576285104275, -1.6508656202000818], [-4.974078211414836, -1.258058633341683, -0.04548698537926161], [-5.454258332708681, -0.30310507905790307, -2.452441772603822], [-6.1535595105033565, 0.6513899495882105, -1.6783822559131625], [-4.823239144821531, 0.0268607394853305, -3.677090739109206], [-6.3950305082002545, -1.5469678375665052, -2.7175632501345497], [-3.5741617863697557, 0.9215323480544237, -1.537966472037064], [-2.6389812609918826, -0.7011341469179394, -2.656638112020167], [-2.0735160409010627, -1.740600756030098, -0.019933065420418592], [-2.5220209528208803, 0.3005938460806614, 0.6115732307470385], [0.0465281587364955, -0.24859461837291646, 0.3246037930626385], [-0.2326570437563611, -0.9037900346867787, -1.7361160911684215], [1.0606513723417526, 1.6513605336852115, -1.7374660375300777], [1.408530802960924, 3.176536224973186, -1.7863201075060056], [2.672413620276367, 3.499806430687654, -2.649509988174067], [3.6650852439076442, 2.6998681533552946, -2.3034648498870838], [2.402646377728109, 3.34562431407816, -3.9297968866850566], [3.025939263676558, 4.752448563851598, -2.429792240096451], [1.6479356871376727, 3.595541709416357, -0.548496662810708], [0.38235591458667256, 3.844326845679148, -2.3112365774402077], [1.8929524931356991, 1.0632336917486571, -0.8851493211874099], [1.2139435881241158, 1.1344810332565651, -2.9535617374347334], [-1.211071099704805, 1.632108284103879, -2.3306200224352915], [-0.7198811239113897, 2.2177451190132307, -0.28051975853154154], [-7.012184668034289, -1.5488900528816725, -3.466341609751753]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'mp2', 'basis': '6-31G', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_R_Addition_MultipleBond_8_5_4', 'label': '5009017845242299296281_R_Addition_MultipleBond_8_5_4', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'guess': 'Mix,Always', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 18 F\n1 28 F\n1 29 F\n2 3 F\n2 16 F\n2 17 F\n3 4 F\n3 14 F\n3 15 F\n4 5 F\n4 12 F\n4 13 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n8 10 F\n8 11 F\n11 30 F\n18 19 F\n18 26 F\n18 27 F\n19 20 F\n19 24 F\n19 25 F\n20 21 F\n20 22 F\n20 23 F\n8 5 4 F\n'}
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
