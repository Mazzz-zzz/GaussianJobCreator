import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_r12_insertion_R_1_5_6'
logfile = '1502984803620600000001_r12_insertion_R_1_5_6.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(-0.9501003540095507), np.float64(-0.1972957561461793), np.float64(-0.27524563161885546)], [np.float64(-1.285673941337389), np.float64(-0.14028930932028008), np.float64(-1.5840991331785415)], [np.float64(-1.420815875440142), np.float64(0.9608602232972487), np.float64(0.22894702084445254)], [np.float64(-1.7425102537906199), np.float64(-1.1524759937293034), np.float64(0.22517447707104826)], [np.float64(0.8128036765422548), np.float64(-0.42473936210214136), np.float64(0.08508089822908144)], [np.float64(0.5279336339740233), np.float64(-0.5432890848512658), np.float64(1.6247460772921687)], [np.float64(1.5087477659288486), np.float64(-1.3099588154383541), np.float64(-0.7717376589380457)], [np.float64(1.3889355171750655), np.float64(1.0676359079406317), np.float64(-0.10607808246635128)], [np.float64(1.1506788317218615), np.float64(1.7395531909945652), np.float64(0.5639780318275157)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r12_insertion_R_1_5_6', 'label': '1502984803620600000001_r12_insertion_R_1_5_6', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 2 F\n1 3 F\n1 4 F\n5 7 F\n5 8 F\n8 9 F\n1 5 F\n5 6 F\n6 1 F\n'}
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
    mol.positions = [[np.float64(-0.9501003540095507), np.float64(-0.1972957561461793), np.float64(-0.27524563161885546)], [np.float64(-1.285673941337389), np.float64(-0.14028930932028008), np.float64(-1.5840991331785415)], [np.float64(-1.420815875440142), np.float64(0.9608602232972487), np.float64(0.22894702084445254)], [np.float64(-1.7425102537906199), np.float64(-1.1524759937293034), np.float64(0.22517447707104826)], [np.float64(0.8128036765422548), np.float64(-0.42473936210214136), np.float64(0.08508089822908144)], [np.float64(0.5279336339740233), np.float64(-0.5432890848512658), np.float64(1.6247460772921687)], [np.float64(1.5087477659288486), np.float64(-1.3099588154383541), np.float64(-0.7717376589380457)], [np.float64(1.3889355171750655), np.float64(1.0676359079406317), np.float64(-0.10607808246635128)], [np.float64(1.1506788317218615), np.float64(1.7395531909945652), np.float64(0.5639780318275157)]]  # reset to the original geometry
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
