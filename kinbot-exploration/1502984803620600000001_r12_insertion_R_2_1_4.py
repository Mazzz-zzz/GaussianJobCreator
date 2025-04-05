import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_r12_insertion_R_2_1_4'
logfile = '1502984803620600000001_r12_insertion_R_2_1_4.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(-0.972602453617316), np.float64(-0.07756844250931584), np.float64(-0.2776661331502726)], [np.float64(-1.3029517761780942), np.float64(-0.19168434206467808), np.float64(-1.6500484848437291)], [np.float64(-1.6640243778788006), np.float64(0.8962734359373395), np.float64(0.20580243481461002)], [np.float64(-1.8398392418957468), np.float64(-1.3197700819932592), np.float64(0.050392992953492234)], [np.float64(0.7497397412290334), np.float64(-0.37369082191453334), np.float64(0.19958613106044357)], [np.float64(0.8236405501359052), np.float64(-0.4511423121956224), np.float64(1.6629237773484975)], [np.float64(1.3695139879257068), np.float64(-1.3131366456818085), np.float64(-0.7176210434033566)], [np.float64(1.3776073842985017), np.float64(1.1075192534239495), np.float64(-0.12056759315314179)], [np.float64(1.4489152638015257), np.float64(1.7231989392754756), np.float64(0.6379609277289868)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r12_insertion_R_2_1_4', 'label': '1502984803620600000001_r12_insertion_R_2_1_4', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 3 F\n1 5 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n2 1 F\n1 4 F\n4 2 F\n'}
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
    mol.positions = [[np.float64(-0.972602453617316), np.float64(-0.07756844250931584), np.float64(-0.2776661331502726)], [np.float64(-1.3029517761780942), np.float64(-0.19168434206467808), np.float64(-1.6500484848437291)], [np.float64(-1.6640243778788006), np.float64(0.8962734359373395), np.float64(0.20580243481461002)], [np.float64(-1.8398392418957468), np.float64(-1.3197700819932592), np.float64(0.050392992953492234)], [np.float64(0.7497397412290334), np.float64(-0.37369082191453334), np.float64(0.19958613106044357)], [np.float64(0.8236405501359052), np.float64(-0.4511423121956224), np.float64(1.6629237773484975)], [np.float64(1.3695139879257068), np.float64(-1.3131366456818085), np.float64(-0.7176210434033566)], [np.float64(1.3776073842985017), np.float64(1.1075192534239495), np.float64(-0.12056759315314179)], [np.float64(1.4489152638015257), np.float64(1.7231989392754756), np.float64(0.6379609277289868)]]  # reset to the original geometry
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
