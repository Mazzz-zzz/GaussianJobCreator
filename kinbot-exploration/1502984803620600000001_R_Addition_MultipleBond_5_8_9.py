import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_R_Addition_MultipleBond_5_8_9'
logfile = '1502984803620600000001_R_Addition_MultipleBond_5_8_9.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(-1.102006555764317), np.float64(-0.19885388144419375), np.float64(-0.3267342210030626)], [np.float64(-1.13446109752887), np.float64(0.049284715519823254), np.float64(-1.63482860363856)], [np.float64(-1.6715755856487926), np.float64(0.7986638622148908), np.float64(0.36263957790634777)], [np.float64(-1.70495761821575), np.float64(-1.3566361231381603), np.float64(-0.053774518444678336)], [np.float64(0.6987941591216608), np.float64(-0.34006948107492035), np.float64(0.21405402482066013)], [np.float64(0.761402459040215), np.float64(-0.3189788186506261), np.float64(1.6680385519493945)], [np.float64(1.3098065023323604), np.float64(-1.393409291682525), np.float64(-0.567063814377968)], [np.float64(1.1810518245375334), np.float64(1.1416578985892525), np.float64(-0.2658873105226633)], [np.float64(1.5332500654404364), np.float64(1.3507115927344933), np.float64(0.7334241647293408)]])

kwargs = {'method': 'mp2', 'basis': '6-31G', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_R_Addition_MultipleBond_5_8_9', 'label': '1502984803620600000001_R_Addition_MultipleBond_5_8_9', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'guess': 'Read,Mix,Always', 'opt': 'ModRedun,Loose,CalcFC', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 6 F\n5 7 F\n5 8 F\n5 8 9 F\n5 8 9 F\n8 9 F\n'}
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
    mol.positions = [[np.float64(-1.102006555764317), np.float64(-0.19885388144419375), np.float64(-0.3267342210030626)], [np.float64(-1.13446109752887), np.float64(0.049284715519823254), np.float64(-1.63482860363856)], [np.float64(-1.6715755856487926), np.float64(0.7986638622148908), np.float64(0.36263957790634777)], [np.float64(-1.70495761821575), np.float64(-1.3566361231381603), np.float64(-0.053774518444678336)], [np.float64(0.6987941591216608), np.float64(-0.34006948107492035), np.float64(0.21405402482066013)], [np.float64(0.761402459040215), np.float64(-0.3189788186506261), np.float64(1.6680385519493945)], [np.float64(1.3098065023323604), np.float64(-1.393409291682525), np.float64(-0.567063814377968)], [np.float64(1.1810518245375334), np.float64(1.1416578985892525), np.float64(-0.2658873105226633)], [np.float64(1.5332500654404364), np.float64(1.3507115927344933), np.float64(0.7334241647293408)]]  # reset to the original geometry
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
