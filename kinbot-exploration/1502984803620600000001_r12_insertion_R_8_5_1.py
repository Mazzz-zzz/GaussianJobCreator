import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_r12_insertion_R_8_5_1'
logfile = '1502984803620600000001_r12_insertion_R_8_5_1.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(-0.9671271029560469), np.float64(-0.17739942330976183), np.float64(-0.2938708630748588)], [np.float64(-1.166076342814378), np.float64(-0.11877837264556713), np.float64(-1.6136507207399218)], [np.float64(-1.4090778516025988), np.float64(0.9995191947091951), np.float64(0.18306456863870701)], [np.float64(-1.8378109141692196), np.float64(-1.0898511154550419), np.float64(0.15445118259920848)], [np.float64(0.8011229221464856), np.float64(-0.48607023623240087), np.float64(0.20920975672814102)], [np.float64(0.8628157689934017), np.float64(-0.5840766740825111), np.float64(1.6537426434968694)], [np.float64(1.4670131966721727), np.float64(-1.3462431065581502), np.float64(-0.7325741033867151)], [np.float64(1.174082954791773), np.float64(1.072068762431076), np.float64(-0.14272134970918973)], [np.float64(1.065056368938411), np.float64(1.7308319711431617), np.float64(0.5731118854477596)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r12_insertion_R_8_5_1', 'label': '1502984803620600000001_r12_insertion_R_8_5_1', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 2 F\n1 3 F\n1 4 F\n5 6 F\n5 7 F\n8 9 F\n8 5 F\n5 1 F\n1 8 F\n'}
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
    mol.positions = [[np.float64(-0.9671271029560469), np.float64(-0.17739942330976183), np.float64(-0.2938708630748588)], [np.float64(-1.166076342814378), np.float64(-0.11877837264556713), np.float64(-1.6136507207399218)], [np.float64(-1.4090778516025988), np.float64(0.9995191947091951), np.float64(0.18306456863870701)], [np.float64(-1.8378109141692196), np.float64(-1.0898511154550419), np.float64(0.15445118259920848)], [np.float64(0.8011229221464856), np.float64(-0.48607023623240087), np.float64(0.20920975672814102)], [np.float64(0.8628157689934017), np.float64(-0.5840766740825111), np.float64(1.6537426434968694)], [np.float64(1.4670131966721727), np.float64(-1.3462431065581502), np.float64(-0.7325741033867151)], [np.float64(1.174082954791773), np.float64(1.072068762431076), np.float64(-0.14272134970918973)], [np.float64(1.065056368938411), np.float64(1.7308319711431617), np.float64(0.5731118854477596)]]  # reset to the original geometry
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
