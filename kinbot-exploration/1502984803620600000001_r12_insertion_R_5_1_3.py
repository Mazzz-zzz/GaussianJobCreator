import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_r12_insertion_R_5_1_3'
logfile = '1502984803620600000001_r12_insertion_R_5_1_3.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(-0.9545893003162812), np.float64(-0.35083622187413827), np.float64(-0.380952880863434)], [np.float64(-1.2651710354676773), np.float64(-0.10450346526584225), np.float64(-1.6201902363320673)], [np.float64(-1.1627870280929382), np.float64(1.0437173813848972), np.float64(0.2740351685431338)], [np.float64(-1.9051100943847743), np.float64(-1.0490694563567178), np.float64(0.17618129771831825)], [np.float64(0.7815365071281809), np.float64(-0.4621421042189572), np.float64(0.19201984155919957)], [np.float64(0.7520743117766977), np.float64(-0.6132844702757296), np.float64(1.6466448175779889)], [np.float64(1.5568571655067116), np.float64(-1.3087964461634232), np.float64(-0.7132035446304542)], [np.float64(1.2517712847671771), np.float64(1.0726000893721714), np.float64(-0.09225358091848347)], [np.float64(0.9354161893291839), np.float64(1.772313693479833), np.float64(0.5084831147188185)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r12_insertion_R_5_1_3', 'label': '1502984803620600000001_r12_insertion_R_5_1_3', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 2 F\n1 4 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n5 1 F\n1 3 F\n3 5 F\n'}
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
    mol.positions = [[np.float64(-0.9545893003162812), np.float64(-0.35083622187413827), np.float64(-0.380952880863434)], [np.float64(-1.2651710354676773), np.float64(-0.10450346526584225), np.float64(-1.6201902363320673)], [np.float64(-1.1627870280929382), np.float64(1.0437173813848972), np.float64(0.2740351685431338)], [np.float64(-1.9051100943847743), np.float64(-1.0490694563567178), np.float64(0.17618129771831825)], [np.float64(0.7815365071281809), np.float64(-0.4621421042189572), np.float64(0.19201984155919957)], [np.float64(0.7520743117766977), np.float64(-0.6132844702757296), np.float64(1.6466448175779889)], [np.float64(1.5568571655067116), np.float64(-1.3087964461634232), np.float64(-0.7132035446304542)], [np.float64(1.2517712847671771), np.float64(1.0726000893721714), np.float64(-0.09225358091848347)], [np.float64(0.9354161893291839), np.float64(1.772313693479833), np.float64(0.5084831147188185)]]  # reset to the original geometry
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
