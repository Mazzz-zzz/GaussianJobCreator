import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_intra_H_migration_suprafacial_6_9'
logfile = '1502984803620600000001_intra_H_migration_suprafacial_6_9.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(2.009341177821025), np.float64(-0.4609842214063651), np.float64(-1.7424174715254614)], [np.float64(3.331424900637044), np.float64(-0.5126274678823691), np.float64(-1.8946216309178927)], [np.float64(1.5744789849417746), np.float64(0.4158154726227546), np.float64(-2.6571872371636203)], [np.float64(1.550589764363197), np.float64(-1.6505212309583936), np.float64(-2.1336516130858088)], [np.float64(1.4554856496314212), np.float64(6.076583246523951e-17), np.float64(-1.228585940635683e-17)], [np.float64(-4.213510553295796e-18), np.float64(2.90752721161415e-18), np.float64(-2.3523962685913924e-17)], [np.float64(2.2583189426830614), np.float64(-0.7264789122146105), np.float64(0.9595269543188124)], [np.float64(1.9885818069143788), np.float64(1.5408568855085618), np.float64(0.0)], [np.float64(1.2210065750421575), np.float64(2.144662637762824), np.float64(-0.0519032763702832)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_intra_H_migration_suprafacial_6_9', 'label': '1502984803620600000001_intra_H_migration_suprafacial_6_9', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n6 5 8 F\n5 8 9 F\n6 5 8 9 F\n'}
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
    mol.positions = [[np.float64(2.009341177821025), np.float64(-0.4609842214063651), np.float64(-1.7424174715254614)], [np.float64(3.331424900637044), np.float64(-0.5126274678823691), np.float64(-1.8946216309178927)], [np.float64(1.5744789849417746), np.float64(0.4158154726227546), np.float64(-2.6571872371636203)], [np.float64(1.550589764363197), np.float64(-1.6505212309583936), np.float64(-2.1336516130858088)], [np.float64(1.4554856496314212), np.float64(6.076583246523951e-17), np.float64(-1.228585940635683e-17)], [np.float64(-4.213510553295796e-18), np.float64(2.90752721161415e-18), np.float64(-2.3523962685913924e-17)], [np.float64(2.2583189426830614), np.float64(-0.7264789122146105), np.float64(0.9595269543188124)], [np.float64(1.9885818069143788), np.float64(1.5408568855085618), np.float64(0.0)], [np.float64(1.2210065750421575), np.float64(2.144662637762824), np.float64(-0.0519032763702832)]]  # reset to the original geometry
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
