import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_r12_insertion_R_6_5_8'
logfile = '1502984803620600000001_r12_insertion_R_6_5_8.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(-1.0359331857934477), np.float64(-0.21041126002294483), np.float64(-0.3193001719739155)], [np.float64(-1.2460621140987473), np.float64(-0.10705703344067977), np.float64(-1.635160547706461)], [np.float64(-1.58020232076836), np.float64(0.8955069823897917), np.float64(0.2111832972411254)], [np.float64(-1.8167573515129085), np.float64(-1.217934192797209), np.float64(0.08583294942429558)], [np.float64(0.7502628057238083), np.float64(-0.4386712356006384), np.float64(0.15654867876980424)], [np.float64(0.852497155738692), np.float64(-0.36498415552487623), np.float64(1.6602595301300347)], [np.float64(1.395812407741768), np.float64(-1.3451528335964305), np.float64(-0.7122085861195855)], [np.float64(1.3150020070934192), np.float64(1.110580844721375), np.float64(-0.12487650003225649)], [np.float64(1.35537959603577), np.float64(1.678122883501445), np.float64(0.668486349982548)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r12_insertion_R_6_5_8', 'label': '1502984803620600000001_r12_insertion_R_6_5_8', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 7 F\n8 9 F\n6 5 F\n5 8 F\n8 6 F\n'}
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
    mol.positions = [[np.float64(-1.0359331857934477), np.float64(-0.21041126002294483), np.float64(-0.3193001719739155)], [np.float64(-1.2460621140987473), np.float64(-0.10705703344067977), np.float64(-1.635160547706461)], [np.float64(-1.58020232076836), np.float64(0.8955069823897917), np.float64(0.2111832972411254)], [np.float64(-1.8167573515129085), np.float64(-1.217934192797209), np.float64(0.08583294942429558)], [np.float64(0.7502628057238083), np.float64(-0.4386712356006384), np.float64(0.15654867876980424)], [np.float64(0.852497155738692), np.float64(-0.36498415552487623), np.float64(1.6602595301300347)], [np.float64(1.395812407741768), np.float64(-1.3451528335964305), np.float64(-0.7122085861195855)], [np.float64(1.3150020070934192), np.float64(1.110580844721375), np.float64(-0.12487650003225649)], [np.float64(1.35537959603577), np.float64(1.678122883501445), np.float64(0.668486349982548)]]  # reset to the original geometry
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
