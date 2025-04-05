import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_r12_insertion_R_6_5_7'
logfile = '1502984803620600000001_r12_insertion_R_6_5_7.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(-1.0716040432265757), np.float64(-0.17978465441412758), np.float64(-0.32107393046037325)], [np.float64(-1.2543811247742973), np.float64(-0.042398290583199244), np.float64(-1.6381211767221402)], [np.float64(-1.686317304155992), np.float64(0.8794076778026345), np.float64(0.23259629738306561)], [np.float64(-1.796617734908725), np.float64(-1.2375567193141452), np.float64(0.03468232819582169)], [np.float64(0.6915922682365644), np.float64(-0.30113903244398527), np.float64(0.17311013599077743)], [np.float64(0.9204589075766294), np.float64(-0.5739528656541667), np.float64(1.639414381728539)], [np.float64(1.474479509841283), np.float64(-1.4131444091967904), np.float64(-0.5385961925077224)], [np.float64(1.2286569662045825), np.float64(1.1343568152087622), np.float64(-0.16016987930917262)], [np.float64(1.4837325705136406), np.float64(1.734210480550979), np.float64(0.5689240402923889)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r12_insertion_R_6_5_7', 'label': '1502984803620600000001_r12_insertion_R_6_5_7', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 8 F\n8 9 F\n6 5 F\n5 7 F\n7 6 F\n'}
Gaussian.command = 'g16 < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy()  # use the Gaussian optimizer
except RuntimeError:
    e = 0.
 
iowait(logfile, 'gauss')
mol.positions = reader_gauss.read_geom(logfile, mol)
if all([ci == 0 for mp in mol.positions for ci in mp]):
    mol.positions = [[np.float64(-1.0716040432265757), np.float64(-0.17978465441412758), np.float64(-0.32107393046037325)], [np.float64(-1.2543811247742973), np.float64(-0.042398290583199244), np.float64(-1.6381211767221402)], [np.float64(-1.686317304155992), np.float64(0.8794076778026345), np.float64(0.23259629738306561)], [np.float64(-1.796617734908725), np.float64(-1.2375567193141452), np.float64(0.03468232819582169)], [np.float64(0.6915922682365644), np.float64(-0.30113903244398527), np.float64(0.17311013599077743)], [np.float64(0.9204589075766294), np.float64(-0.5739528656541667), np.float64(1.639414381728539)], [np.float64(1.474479509841283), np.float64(-1.4131444091967904), np.float64(-0.5385961925077224)], [np.float64(1.2286569662045825), np.float64(1.1343568152087622), np.float64(-0.16016987930917262)], [np.float64(1.4837325705136406), np.float64(1.734210480550979), np.float64(0.5689240402923889)]]  # reset to the original geometry
db.write(mol, name=label, data={'energy': e, 'status': 'normal'})

#for tr in range(ntrial):  # DELETED CURLY BRACKET
#    try:
#        success = True
#        e = mol.get_potential_energy() # use the Gaussian optimizer (task optimize)
#        iowait(logfile, 'gauss')
#        mol.positions = reader_gauss.read_geom(logfile, mol)
#        db.write(mol, name=label, data={'energy': e,'status': 'normal'})
#        break
#    except RuntimeError: 
#        success = False
#        
#if not success:
#    if not bimol:
#        try:
#            mol.positions = reader_gauss.read_geom(logfile, mol)
#            del kwargs['opt']  # this is when we give up optimization!!
#            calc = Gaussian(**kwargs)
#            e = mol.get_potential_energy() 
#            iowait(logfile, 'gauss')
#            mol.positions = reader_gauss.read_geom(logfile, mol)
#            db.write(mol, name=label, data={'energy': e,'status': 'normal'})
#        except: 
#            db.write(mol, name = label, data = {'status': 'error'})
#    else:
#        try:
#            mol.positions = reader_gauss.read_geom(logfile, mol)
#            db.write(mol, name=label, data={'energy': e,'status': 'normal'})
#        except: 
#            db.write(mol, name = label, data = {'status': 'error'})

with open(logfile,'a') as f:
    f.write('done\n')
