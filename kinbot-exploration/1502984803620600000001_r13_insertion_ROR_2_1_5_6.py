import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_r13_insertion_ROR_2_1_5_6'
logfile = '1502984803620600000001_r13_insertion_ROR_2_1_5_6.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(1.331817288232136), np.float64(-1.1646806534043306e-25), np.float64(-7.085689672884592e-26)], [np.float64(-2.174527528367906e-23), np.float64(-1.1724250888101834e-23), np.float64(-6.20328958521337e-24)], [np.float64(1.6862625258972894), np.float64(-0.7162514848332132), np.float64(1.075191102066377)], [np.float64(1.701627498305902), np.float64(-0.7269608453514205), np.float64(-1.055122356601425)], [np.float64(2.089742862227282), np.float64(1.7265056840491804), np.float64(4.458051554712986e-26)], [np.float64(3.211149929880731), np.float64(1.7313891932669274), np.float64(0.9278268573929153)], [np.float64(2.110723304839141), np.float64(2.2242772840048595), np.float64(-1.3582311148366368)], [np.float64(0.9242423321676769), np.float64(2.5283002739415594), np.float64(0.8106562878829063)], [np.float64(1.0817030548873088), np.float64(2.7092417363828343), np.float64(1.75876786993442)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r13_insertion_ROR_2_1_5_6', 'label': '1502984803620600000001_r13_insertion_ROR_2_1_5_6', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 2 F\n1 3 F\n1 4 F\n1 5 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n2 1 5 6 F\n'}
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
    mol.positions = [[np.float64(1.331817288232136), np.float64(-1.1646806534043306e-25), np.float64(-7.085689672884592e-26)], [np.float64(-2.174527528367906e-23), np.float64(-1.1724250888101834e-23), np.float64(-6.20328958521337e-24)], [np.float64(1.6862625258972894), np.float64(-0.7162514848332132), np.float64(1.075191102066377)], [np.float64(1.701627498305902), np.float64(-0.7269608453514205), np.float64(-1.055122356601425)], [np.float64(2.089742862227282), np.float64(1.7265056840491804), np.float64(4.458051554712986e-26)], [np.float64(3.211149929880731), np.float64(1.7313891932669274), np.float64(0.9278268573929153)], [np.float64(2.110723304839141), np.float64(2.2242772840048595), np.float64(-1.3582311148366368)], [np.float64(0.9242423321676769), np.float64(2.5283002739415594), np.float64(0.8106562878829063)], [np.float64(1.0817030548873088), np.float64(2.7092417363828343), np.float64(1.75876786993442)]]  # reset to the original geometry
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
