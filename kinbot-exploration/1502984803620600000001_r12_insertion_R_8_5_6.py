import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_r12_insertion_R_8_5_6'
logfile = '1502984803620600000001_r12_insertion_R_8_5_6.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(-1.0369106942773996), np.float64(-0.20822841718937576), np.float64(-0.3220185277483905)], [np.float64(-1.254344882573616), np.float64(-0.10740180071696728), np.float64(-1.6385032033530407)], [np.float64(-1.5830544697280085), np.float64(0.8967007866275899), np.float64(0.20917433233687577)], [np.float64(-1.8112136377552324), np.float64(-1.2197109123997052), np.float64(0.0883247514988148)], [np.float64(0.7550298817964999), np.float64(-0.4264095749864349), np.float64(0.1422631732818485)], [np.float64(0.8487869966760541), np.float64(-0.3682297248488768), np.float64(1.6798562432268695)], [np.float64(1.3980988396586502), np.float64(-1.3400724659146506), np.float64(-0.7170345494833327)], [np.float64(1.314574376961781), np.float64(1.0937622237882345), np.float64(-0.11557942534859031)], [np.float64(1.3590345901176126), np.float64(1.6795898844352188), np.float64(0.6642802087656763)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r12_insertion_R_8_5_6', 'label': '1502984803620600000001_r12_insertion_R_8_5_6', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 7 F\n8 9 F\n8 5 F\n5 6 F\n6 8 F\n'}
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
    mol.positions = [[np.float64(-1.0369106942773996), np.float64(-0.20822841718937576), np.float64(-0.3220185277483905)], [np.float64(-1.254344882573616), np.float64(-0.10740180071696728), np.float64(-1.6385032033530407)], [np.float64(-1.5830544697280085), np.float64(0.8967007866275899), np.float64(0.20917433233687577)], [np.float64(-1.8112136377552324), np.float64(-1.2197109123997052), np.float64(0.0883247514988148)], [np.float64(0.7550298817964999), np.float64(-0.4264095749864349), np.float64(0.1422631732818485)], [np.float64(0.8487869966760541), np.float64(-0.3682297248488768), np.float64(1.6798562432268695)], [np.float64(1.3980988396586502), np.float64(-1.3400724659146506), np.float64(-0.7170345494833327)], [np.float64(1.314574376961781), np.float64(1.0937622237882345), np.float64(-0.11557942534859031)], [np.float64(1.3590345901176126), np.float64(1.6795898844352188), np.float64(0.6642802087656763)]]  # reset to the original geometry
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
