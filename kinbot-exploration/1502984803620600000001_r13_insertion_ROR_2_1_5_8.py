import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_r13_insertion_ROR_2_1_5_8'
logfile = '1502984803620600000001_r13_insertion_ROR_2_1_5_8.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(1.3318173925857855), np.float64(-5.81001135870373e-17), np.float64(-3.620653810150714e-17)], [np.float64(-2.7686984322489385e-17), np.float64(-8.879868627437933e-17), np.float64(3.3751627335847897e-17)], [np.float64(1.6882450375562075), np.float64(-0.7138250839366063), np.float64(1.076150222852152)], [np.float64(1.7001324289749837), np.float64(-0.7273266573883941), np.float64(-1.0553948111099702)], [np.float64(2.0940306346850495), np.float64(1.7246182462033528), np.float64(-2.9277665780460695e-17)], [np.float64(3.473217087847918), np.float64(1.6074799977368375), np.float64(0.4500665668578651)], [np.float64(1.647692260450913), np.float64(2.4306964837385556), np.float64(-1.1812034875592814)], [np.float64(1.3122203000479384), np.float64(2.356009222391638), np.float64(1.2839587828815269)], [np.float64(1.7841976766804963), np.float64(2.3593995229351146), np.float64(2.1405071755023255)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r13_insertion_ROR_2_1_5_8', 'label': '1502984803620600000001_r13_insertion_ROR_2_1_5_8', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n2 1 5 8 F\n'}
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
    mol.positions = [[np.float64(1.3318173925857855), np.float64(-5.81001135870373e-17), np.float64(-3.620653810150714e-17)], [np.float64(-2.7686984322489385e-17), np.float64(-8.879868627437933e-17), np.float64(3.3751627335847897e-17)], [np.float64(1.6882450375562075), np.float64(-0.7138250839366063), np.float64(1.076150222852152)], [np.float64(1.7001324289749837), np.float64(-0.7273266573883941), np.float64(-1.0553948111099702)], [np.float64(2.0940306346850495), np.float64(1.7246182462033528), np.float64(-2.9277665780460695e-17)], [np.float64(3.473217087847918), np.float64(1.6074799977368375), np.float64(0.4500665668578651)], [np.float64(1.647692260450913), np.float64(2.4306964837385556), np.float64(-1.1812034875592814)], [np.float64(1.3122203000479384), np.float64(2.356009222391638), np.float64(1.2839587828815269)], [np.float64(1.7841976766804963), np.float64(2.3593995229351146), np.float64(2.1405071755023255)]]  # reset to the original geometry
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
