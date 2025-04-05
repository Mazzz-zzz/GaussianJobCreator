import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_r12_insertion_R_2_1_3'
logfile = '1502984803620600000001_r12_insertion_R_2_1_3.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(-0.9066953312934298), np.float64(-0.36795727445069565), np.float64(-0.24697098783416885)], [np.float64(-1.2140482784869717), np.float64(0.03359369136822543), np.float64(-1.5999485635698307)], [np.float64(-1.378946112707795), np.float64(1.1322113544401091), np.float64(0.1670626778752646)], [np.float64(-1.8595165636799058), np.float64(-1.0814473016924608), np.float64(0.1696335926973728)], [np.float64(0.8198703396535482), np.float64(-0.4987435754268192), np.float64(0.21808683618491087)], [np.float64(0.8587311474165776), np.float64(-0.667894413097862), np.float64(1.6798118638757313)], [np.float64(1.5323880935805039), np.float64(-1.3572179999465641), np.float64(-0.7373632455494591)], [np.float64(1.2835325622119207), np.float64(1.041214361187052), np.float64(-0.08123290718952424)], [np.float64(0.8546841432515565), np.float64(1.7662411371044215), np.float64(0.42168672499455967)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r12_insertion_R_2_1_3', 'label': '1502984803620600000001_r12_insertion_R_2_1_3', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 4 F\n1 5 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n2 1 F\n1 3 F\n3 2 F\n'}
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
    mol.positions = [[np.float64(-0.9066953312934298), np.float64(-0.36795727445069565), np.float64(-0.24697098783416885)], [np.float64(-1.2140482784869717), np.float64(0.03359369136822543), np.float64(-1.5999485635698307)], [np.float64(-1.378946112707795), np.float64(1.1322113544401091), np.float64(0.1670626778752646)], [np.float64(-1.8595165636799058), np.float64(-1.0814473016924608), np.float64(0.1696335926973728)], [np.float64(0.8198703396535482), np.float64(-0.4987435754268192), np.float64(0.21808683618491087)], [np.float64(0.8587311474165776), np.float64(-0.667894413097862), np.float64(1.6798118638757313)], [np.float64(1.5323880935805039), np.float64(-1.3572179999465641), np.float64(-0.7373632455494591)], [np.float64(1.2835325622119207), np.float64(1.041214361187052), np.float64(-0.08123290718952424)], [np.float64(0.8546841432515565), np.float64(1.7662411371044215), np.float64(0.42168672499455967)]]  # reset to the original geometry
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
