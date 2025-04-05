import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_intra_H_migration_6_9'
logfile = '1502984803620600000001_intra_H_migration_6_9.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(2.038499330070244), np.float64(-0.447145004885665), np.float64(-1.7364997225168122)], [np.float64(3.3620783598171835), np.float64(-0.43077420942931105), np.float64(-1.8834937254022237)], [np.float64(1.5615104634943102), np.float64(0.39213553679023294), np.float64(-2.665357144590199)], [np.float64(1.644711535476605), np.float64(-1.664722344162494), np.float64(-2.111934060596901)], [np.float64(1.4554852311851192), np.float64(1.4198372181225923e-18), np.float64(-1.7608512290924065e-20)], [np.float64(6.171323226987143e-17), np.float64(1.230368288330406e-16), np.float64(-3.779981917718811e-18)], [np.float64(2.297346028244733), np.float64(-0.6733423498040921), np.float64(0.964825865303919)], [np.float64(1.853135589669854), np.float64(1.58123569341922), np.float64(-3.0966766352439785e-19)], [np.float64(1.1236014399112348), np.float64(2.2292974232904363), np.float64(-0.0651509708776127)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_intra_H_migration_6_9', 'label': '1502984803620600000001_intra_H_migration_6_9', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n6 5 8 F\n5 8 9 F\n6 5 8 9 F\n'}
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
    mol.positions = [[np.float64(2.038499330070244), np.float64(-0.447145004885665), np.float64(-1.7364997225168122)], [np.float64(3.3620783598171835), np.float64(-0.43077420942931105), np.float64(-1.8834937254022237)], [np.float64(1.5615104634943102), np.float64(0.39213553679023294), np.float64(-2.665357144590199)], [np.float64(1.644711535476605), np.float64(-1.664722344162494), np.float64(-2.111934060596901)], [np.float64(1.4554852311851192), np.float64(1.4198372181225923e-18), np.float64(-1.7608512290924065e-20)], [np.float64(6.171323226987143e-17), np.float64(1.230368288330406e-16), np.float64(-3.779981917718811e-18)], [np.float64(2.297346028244733), np.float64(-0.6733423498040921), np.float64(0.964825865303919)], [np.float64(1.853135589669854), np.float64(1.58123569341922), np.float64(-3.0966766352439785e-19)], [np.float64(1.1236014399112348), np.float64(2.2292974232904363), np.float64(-0.0651509708776127)]]  # reset to the original geometry
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
