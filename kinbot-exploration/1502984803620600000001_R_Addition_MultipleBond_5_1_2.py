import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_R_Addition_MultipleBond_5_1_2'
logfile = '1502984803620600000001_R_Addition_MultipleBond_5_1_2.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(-1.0755808895295689), np.float64(-0.20929654817853136), np.float64(-0.2568514169751835)], [np.float64(-0.6864921521536994), np.float64(-0.039310204163839325), np.float64(-1.53361093506564)], [np.float64(-1.7227865490079788), np.float64(0.8517138344506456), np.float64(0.2390137355429403)], [np.float64(-1.777779675790621), np.float64(-1.3286595410369084), np.float64(-0.08920328574185886)], [np.float64(0.727033464303939), np.float64(-0.35754044987826844), np.float64(0.27241797720325545)], [np.float64(0.8756007373498693), np.float64(-0.4617568527314058), np.float64(1.7161127300501873)], [np.float64(1.3519527341888427), np.float64(-1.3184387999512903), np.float64(-0.6100235128697794)], [np.float64(1.09063937285341), np.float64(1.1716850809866057), np.float64(-0.1604346078062861)], [np.float64(1.6726856002552464), np.float64(1.6505373797650418), np.float64(0.4625984243812675)]])

kwargs = {'method': 'mp2', 'basis': '6-31G', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_R_Addition_MultipleBond_5_1_2', 'label': '1502984803620600000001_R_Addition_MultipleBond_5_1_2', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'guess': 'Read,Mix,Always', 'opt': 'ModRedun,Loose,CalcFC', 'addsec': '1 3 F\n1 4 F\n1 5 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n5 1 2 F\n5 1 2 F\n1 2 F\n'}
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
    mol.positions = [[np.float64(-1.0755808895295689), np.float64(-0.20929654817853136), np.float64(-0.2568514169751835)], [np.float64(-0.6864921521536994), np.float64(-0.039310204163839325), np.float64(-1.53361093506564)], [np.float64(-1.7227865490079788), np.float64(0.8517138344506456), np.float64(0.2390137355429403)], [np.float64(-1.777779675790621), np.float64(-1.3286595410369084), np.float64(-0.08920328574185886)], [np.float64(0.727033464303939), np.float64(-0.35754044987826844), np.float64(0.27241797720325545)], [np.float64(0.8756007373498693), np.float64(-0.4617568527314058), np.float64(1.7161127300501873)], [np.float64(1.3519527341888427), np.float64(-1.3184387999512903), np.float64(-0.6100235128697794)], [np.float64(1.09063937285341), np.float64(1.1716850809866057), np.float64(-0.1604346078062861)], [np.float64(1.6726856002552464), np.float64(1.6505373797650418), np.float64(0.4625984243812675)]]  # reset to the original geometry
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
