import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_r13_insertion_ROR_9_8_5_6'
logfile = '1502984803620600000001_r13_insertion_ROR_9_8_5_6.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(2.269743629230589), np.float64(1.4990411660374456), np.float64(1.8043520069235799)], [np.float64(3.2414812779003266), np.float64(0.6376160277539176), np.float64(2.0999997788651465)], [np.float64(1.266436358803968), np.float64(1.2140824395708827), np.float64(2.6451006669276484)], [np.float64(2.72133291542358), np.float64(2.6987321887205287), np.float64(2.172221552756503)], [np.float64(1.7246468705386495), np.float64(1.4494549721585843), np.float64(1.7688027489748718e-18)], [np.float64(0.6363591141920893), np.float64(2.4015587356184223), np.float64(-0.1660282596270403)], [np.float64(2.8970314739286644), np.float64(1.3509427619936025), np.float64(-0.8419116651249388)], [np.float64(0.9779824015395164), np.float64(-1.4358830337056634e-18), np.float64(-5.930155602684091e-19)], [np.float64(5.572954724541562e-17), np.float64(2.1031800760676714e-16), np.float64(-1.454002862944191e-17)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r13_insertion_ROR_9_8_5_6', 'label': '1502984803620600000001_r13_insertion_ROR_9_8_5_6', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n9 8 5 6 F\n'}
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
    mol.positions = [[np.float64(2.269743629230589), np.float64(1.4990411660374456), np.float64(1.8043520069235799)], [np.float64(3.2414812779003266), np.float64(0.6376160277539176), np.float64(2.0999997788651465)], [np.float64(1.266436358803968), np.float64(1.2140824395708827), np.float64(2.6451006669276484)], [np.float64(2.72133291542358), np.float64(2.6987321887205287), np.float64(2.172221552756503)], [np.float64(1.7246468705386495), np.float64(1.4494549721585843), np.float64(1.7688027489748718e-18)], [np.float64(0.6363591141920893), np.float64(2.4015587356184223), np.float64(-0.1660282596270403)], [np.float64(2.8970314739286644), np.float64(1.3509427619936025), np.float64(-0.8419116651249388)], [np.float64(0.9779824015395164), np.float64(-1.4358830337056634e-18), np.float64(-5.930155602684091e-19)], [np.float64(5.572954724541562e-17), np.float64(2.1031800760676714e-16), np.float64(-1.454002862944191e-17)]]  # reset to the original geometry
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
