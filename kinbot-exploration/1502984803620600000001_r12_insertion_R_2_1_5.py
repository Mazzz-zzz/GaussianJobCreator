import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/kinbot.db')
label = '1502984803620600000001_r12_insertion_R_2_1_5'
logfile = '1502984803620600000001_r12_insertion_R_2_1_5.log'

scan = 0
bimol = 0
mol = Atoms(symbols=[np.str_('C'), np.str_('F'), np.str_('F'), np.str_('F'), np.str_('S'), np.str_('O'), np.str_('O'), np.str_('O'), np.str_('H')], positions=[[np.float64(-1.1345783529783418), np.float64(-0.2128154756606224), np.float64(-0.16462441559355448)], [np.float64(-0.8932519962739979), np.float64(0.0115165557422985), np.float64(-1.542088300769852)], [np.float64(-1.7297979755304869), np.float64(0.8424680549984556), np.float64(0.2917914714101018)], [np.float64(-1.8886735041001927), np.float64(-1.2515666156682748), np.float64(-0.06877266676219534)], [np.float64(0.7805436791576351), np.float64(-0.4097745917912974), np.float64(0.22426572282213414)], [np.float64(0.9300737267907326), np.float64(-0.42166593201672126), np.float64(1.6656457309423527)], [np.float64(1.3103011002305351), np.float64(-1.374817679315529), np.float64(-0.6940039920683378)], [np.float64(1.252224248731643), np.float64(1.073258183720408), np.float64(-0.21141762147434404)], [np.float64(1.3631580739724711), np.float64(1.7433984999912804), np.float64(0.4899660714936957)]])

kwargs = {'method': 'am1', 'basis': '', 'nprocshared': 8, 'mem': '700MW', 'chk': '1502984803620600000001_r12_insertion_R_2_1_5', 'label': '1502984803620600000001_r12_insertion_R_2_1_5', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'opt': 'ModRedun,Loose,CalcFC', 'guess': 'Read', 'addsec': '1 3 F\n1 4 F\n5 6 F\n5 7 F\n5 8 F\n8 9 F\n2 1 F\n1 5 F\n5 2 F\n'}
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
    mol.positions = [[np.float64(-1.1345783529783418), np.float64(-0.2128154756606224), np.float64(-0.16462441559355448)], [np.float64(-0.8932519962739979), np.float64(0.0115165557422985), np.float64(-1.542088300769852)], [np.float64(-1.7297979755304869), np.float64(0.8424680549984556), np.float64(0.2917914714101018)], [np.float64(-1.8886735041001927), np.float64(-1.2515666156682748), np.float64(-0.06877266676219534)], [np.float64(0.7805436791576351), np.float64(-0.4097745917912974), np.float64(0.22426572282213414)], [np.float64(0.9300737267907326), np.float64(-0.42166593201672126), np.float64(1.6656457309423527)], [np.float64(1.3103011002305351), np.float64(-1.374817679315529), np.float64(-0.6940039920683378)], [np.float64(1.252224248731643), np.float64(1.073258183720408), np.float64(-0.21141762147434404)], [np.float64(1.3631580739724711), np.float64(1.7433984999912804), np.float64(0.4899660714936957)]]  # reset to the original geometry
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
