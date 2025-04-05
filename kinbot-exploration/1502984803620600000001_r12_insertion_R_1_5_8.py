import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_r12_insertion_R_1_5_8'
logfile = '1502984803620600000001_r12_insertion_R_1_5_8.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(-0.9320202043597786), np.float64(-0.1762444132041466), np.float64(-0.2818627654560487)], [np.float64(-1.1308942788562373), np.float64(-0.07493963149786906), np.float64(-1.6054638800014351)], [np.float64(-1.3865459607125619), np.float64(0.9916885069769944), np.float64(0.2215956352294351)], [np.float64(-1.837114461413759), np.float64(-1.089547889703256), np.float64(0.12290702076739574)], [np.float64(0.778969649638584), np.float64(-0.5146606457289055), np.float64(0.21562229557044682)], [np.float64(0.8755562603834754), np.float64(-0.5738119301483847), np.float64(1.6587747415486034)], [np.float64(1.4583628107596958), np.float64(-1.3758984022566543), np.float64(-0.7116868850484256)], [np.float64(1.128273536845004), np.float64(1.0871165067488746), np.float64(-0.18213404262424512)], [np.float64(1.035414647430454), np.float64(1.7262998983975166), np.float64(0.5530118803418383)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r12_insertion_R_1_5_8', 'label': '1502984803620600000001_r12_insertion_R_1_5_8', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 2 F\n1 3 F\n1 4 F\n5 6 F\n5 7 F\n8 9 F\n1 5 F\n5 8 F\n8 1 F\n'}
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
    mol.positions = [[np.float64(-0.9320202043597786), np.float64(-0.1762444132041466), np.float64(-0.2818627654560487)], [np.float64(-1.1308942788562373), np.float64(-0.07493963149786906), np.float64(-1.6054638800014351)], [np.float64(-1.3865459607125619), np.float64(0.9916885069769944), np.float64(0.2215956352294351)], [np.float64(-1.837114461413759), np.float64(-1.089547889703256), np.float64(0.12290702076739574)], [np.float64(0.778969649638584), np.float64(-0.5146606457289055), np.float64(0.21562229557044682)], [np.float64(0.8755562603834754), np.float64(-0.5738119301483847), np.float64(1.6587747415486034)], [np.float64(1.4583628107596958), np.float64(-1.3758984022566543), np.float64(-0.7116868850484256)], [np.float64(1.128273536845004), np.float64(1.0871165067488746), np.float64(-0.18213404262424512)], [np.float64(1.035414647430454), np.float64(1.7262998983975166), np.float64(0.5530118803418383)]]  # reset to the original geometry
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
