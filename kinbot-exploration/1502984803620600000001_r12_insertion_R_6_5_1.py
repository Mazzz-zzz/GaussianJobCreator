import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_r12_insertion_R_6_5_1'
logfile = '1502984803620600000001_r12_insertion_R_6_5_1.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(-0.9793661526927114), np.float64(-0.19214938797321895), np.float64(-0.28932064808604413)], [np.float64(-1.277510850558725), np.float64(-0.15178075795137963), np.float64(-1.5951582063971221)], [np.float64(-1.4418720419906408), np.float64(0.9674562612532074), np.float64(0.206386014374987)], [np.float64(-1.7684354654465049), np.float64(-1.1411412856680285), np.float64(0.21875465665752217)], [np.float64(0.8230556727698369), np.float64(-0.430878584299265), np.float64(0.12206612365999646)], [np.float64(0.6026221859462672), np.float64(-0.5579984977889507), np.float64(1.6095089802376998)], [np.float64(1.5003700398341706), np.float64(-1.3124938764131886), np.float64(-0.7619504856934273)], [np.float64(1.3766918858075181), np.float64(1.0699718298141387), np.float64(-0.09480056676700264)], [np.float64(1.1544447265596294), np.float64(1.749013298540398), np.float64(0.57527713193711)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r12_insertion_R_6_5_1', 'label': '1502984803620600000001_r12_insertion_R_6_5_1', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 2 F\n1 3 F\n1 4 F\n5 7 F\n5 8 F\n8 9 F\n6 5 F\n5 1 F\n1 6 F\n'}
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
    mol.positions = [[np.float64(-0.9793661526927114), np.float64(-0.19214938797321895), np.float64(-0.28932064808604413)], [np.float64(-1.277510850558725), np.float64(-0.15178075795137963), np.float64(-1.5951582063971221)], [np.float64(-1.4418720419906408), np.float64(0.9674562612532074), np.float64(0.206386014374987)], [np.float64(-1.7684354654465049), np.float64(-1.1411412856680285), np.float64(0.21875465665752217)], [np.float64(0.8230556727698369), np.float64(-0.430878584299265), np.float64(0.12206612365999646)], [np.float64(0.6026221859462672), np.float64(-0.5579984977889507), np.float64(1.6095089802376998)], [np.float64(1.5003700398341706), np.float64(-1.3124938764131886), np.float64(-0.7619504856934273)], [np.float64(1.3766918858075181), np.float64(1.0699718298141387), np.float64(-0.09480056676700264)], [np.float64(1.1544447265596294), np.float64(1.749013298540398), np.float64(0.57527713193711)]]  # reset to the original geometry
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
