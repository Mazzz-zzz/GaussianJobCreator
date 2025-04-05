import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_r12_insertion_R_9_8_5'
logfile = '1502984803620600000001_r12_insertion_R_9_8_5.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(-1.0450987322757497), np.float64(-0.19410435983481708), np.float64(-0.3224700437434299)], [np.float64(-1.2537748645900135), np.float64(-0.11839452459681732), np.float64(-1.635253793644243)], [np.float64(-1.6080053707072284), np.float64(0.9031644533414231), np.float64(0.1957719018230253)], [np.float64(-1.7724621049163465), np.float64(-1.2196179599910784), np.float64(0.10985132207136944)], [np.float64(0.789113389241341), np.float64(-0.33230758133398913), np.float64(0.1465900070493189)], [np.float64(0.8386179132508352), np.float64(-0.36602103081584847), np.float64(1.6222945795026962)], [np.float64(1.3841296055327477), np.float64(-1.3678749687538048), np.float64(-0.7063571146319292)], [np.float64(1.2963513649128697), np.float64(1.0323576052399468), np.float64(-0.19471097089148284)], [np.float64(1.3611287975084283), np.float64(1.6628003669003333), np.float64(0.7750491134466501)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r12_insertion_R_9_8_5', 'label': '1502984803620600000001_r12_insertion_R_9_8_5', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 6 F\n5 7 F\n9 8 F\n8 5 F\n5 9 F\n'}
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
    mol.positions = [[np.float64(-1.0450987322757497), np.float64(-0.19410435983481708), np.float64(-0.3224700437434299)], [np.float64(-1.2537748645900135), np.float64(-0.11839452459681732), np.float64(-1.635253793644243)], [np.float64(-1.6080053707072284), np.float64(0.9031644533414231), np.float64(0.1957719018230253)], [np.float64(-1.7724621049163465), np.float64(-1.2196179599910784), np.float64(0.10985132207136944)], [np.float64(0.789113389241341), np.float64(-0.33230758133398913), np.float64(0.1465900070493189)], [np.float64(0.8386179132508352), np.float64(-0.36602103081584847), np.float64(1.6222945795026962)], [np.float64(1.3841296055327477), np.float64(-1.3678749687538048), np.float64(-0.7063571146319292)], [np.float64(1.2963513649128697), np.float64(1.0323576052399468), np.float64(-0.19471097089148284)], [np.float64(1.3611287975084283), np.float64(1.6628003669003333), np.float64(0.7750491134466501)]]  # reset to the original geometry
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
