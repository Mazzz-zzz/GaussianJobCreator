import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_r12_insertion_R_5_8_9'
logfile = '1502984803620600000001_r12_insertion_R_5_8_9.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(-1.0214093889674827), np.float64(-0.19475968030129467), np.float64(-0.30772547579582865)], [np.float64(-1.2133499345651884), np.float64(-0.10322248160367076), np.float64(-1.6220503548016887)], [np.float64(-1.551471039742581), np.float64(0.9249278499518507), np.float64(0.21042054257159015)], [np.float64(-1.8125496763909037), np.float64(-1.188630138225831), np.float64(0.10359222441975538)], [np.float64(0.7654094843380119), np.float64(-0.4257299723209253), np.float64(0.19051723211415197)], [np.float64(0.8636095736364713), np.float64(-0.3649527625713504), np.float64(1.6229075899909156)], [np.float64(1.4149619234922572), np.float64(-1.3441757385501292), np.float64(-0.6864895755420513)], [np.float64(1.3092442612488586), np.float64(1.2055110056266929), np.float64(-0.24110239822391755)], [np.float64(1.2355537988177319), np.float64(1.4910329197729206), np.float64(0.720694215533812)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r12_insertion_R_5_8_9', 'label': '1502984803620600000001_r12_insertion_R_5_8_9', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n9 5 F\n'}
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
    mol.positions = [[np.float64(-1.0214093889674827), np.float64(-0.19475968030129467), np.float64(-0.30772547579582865)], [np.float64(-1.2133499345651884), np.float64(-0.10322248160367076), np.float64(-1.6220503548016887)], [np.float64(-1.551471039742581), np.float64(0.9249278499518507), np.float64(0.21042054257159015)], [np.float64(-1.8125496763909037), np.float64(-1.188630138225831), np.float64(0.10359222441975538)], [np.float64(0.7654094843380119), np.float64(-0.4257299723209253), np.float64(0.19051723211415197)], [np.float64(0.8636095736364713), np.float64(-0.3649527625713504), np.float64(1.6229075899909156)], [np.float64(1.4149619234922572), np.float64(-1.3441757385501292), np.float64(-0.6864895755420513)], [np.float64(1.3092442612488586), np.float64(1.2055110056266929), np.float64(-0.24110239822391755)], [np.float64(1.2355537988177319), np.float64(1.4910329197729206), np.float64(0.720694215533812)]]  # reset to the original geometry
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
