import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = '5009017845242299296281_R_Addition_MultipleBond_8_11_30'
logfile = '5009017845242299296281_R_Addition_MultipleBond_8_11_30.log'

atom = ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H']
geom = [[0.034248836183416974, 0.08850266492249827, -0.00098467502159751], [-0.6839690693511622, -1.284545641566585, 0.2007396650142067], [-2.261003210106903, -1.2460606310044424, 0.18626688301617714], [-2.9757119322930063, -0.7468147297548774, -1.120878582843062], [-4.432331448321425, -1.3165015821950252, -1.2290033816050097], [-4.398524587039332, -2.5505056191865254, -1.7204589018871397], [-5.0391200082191245, -1.269454252290657, -0.053657489582170516], [-5.487429497791386, -0.3013310514944793, -2.404795730554491], [-5.869641363802847, 0.8890620714621034, -1.7427416291256654], [-4.8372866375900365, -0.3131481135611835, -3.670598053624814], [-6.7045449655507126, -1.290769997709783, -2.535712895419884], [-3.066070513180156, 0.5793268302970047, -1.0922826034020816], [-2.2955172852434105, -1.1503357626648056, -2.187070609420142], [-2.6386903638997037, -2.513260716744795, 0.38115549580013164], [-2.662168715308759, -0.4868003703609385, 1.1986714826391944], [-0.3611074883625628, -1.750356123187119, 1.4052012336081823], [-0.26982924069476233, -2.128581988458272, -0.7411188697750919], [1.5902534135528443, -0.055760089718434595, -0.17768987928746735], [2.341813462661274, 1.2947075913788728, 0.08771469145064983], [3.7931683740386513, 1.3234020332119365, -0.49785783063583894], [4.436578002988003, 0.22598973442816966, -0.14190951882582517], [3.7525058050369147, 1.4106736961084303, -1.8114998551319892], [4.421663592686021, 2.37628071384169, -0.012062249867792378], [2.429719425573684, 1.4768121191460968, 1.4005543629094614], [1.659587682369306, 2.2897500580446333, -0.47585629684951897], [2.0493218158657203, -0.96335261597485, 0.6785321023811502], [1.829804387863204, -0.443152553906519, -1.4284735404543918], [-0.44964392133810543, 0.6560752733555626, -1.103075900548324], [-0.21659948780153288, 0.8484498314195642, 1.0625176497266184], [-6.357519139508325, -1.5433898445200422, -3.852718085262044]]
mol = Atoms(symbols=atom, positions=geom)

kwargs = {'method': 'mp2', 'basis': '6-31G', 'nprocshared': 8, 'mem': '700MW', 'chk': '5009017845242299296281_R_Addition_MultipleBond_8_11_30', 'label': '5009017845242299296281_R_Addition_MultipleBond_8_11_30', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'guess': 'Mix,Always', 'opt': 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999', 'addsec': '1 2 F\n1 18 F\n1 28 F\n1 29 F\n2 3 F\n2 16 F\n2 17 F\n3 4 F\n3 14 F\n3 15 F\n4 5 F\n4 12 F\n4 13 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n8 10 F\n8 11 F\n11 30 F\n18 19 F\n18 26 F\n18 27 F\n19 20 F\n19 24 F\n19 25 F\n20 21 F\n20 22 F\n20 23 F\n8 11 30 F\n'}
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
