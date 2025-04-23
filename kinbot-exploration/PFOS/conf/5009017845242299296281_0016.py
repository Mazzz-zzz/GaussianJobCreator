import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0016'
logfile = 'conf/5009017845242299296281_0016.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586384, -1.3935598728845964, 0.0866492574076532], [-0.34660204151390517, -2.4204394252486656, -1.0595513872112563], [1.1624535490467045, -2.7999021924931125, -1.2603960465750406], [1.322955526797157, -4.171696370783136, -2.003366440702424], [1.1329590874195947, -5.166849009408313, -1.1547943919415717], [0.4456504601938104, -4.264656044096986, -2.998652510032487], [3.0339402189610625, -4.3645246490660625, -2.749698785834422], [3.9877413956025296, -3.905181555060063, -1.8001208557853121], [3.0744314826683046, -5.641324489848281, -3.360844755133002], [2.904675852559841, -3.2914695967279375, -3.895293045613823], [1.7500514572469128, -1.8550381832545377, -1.9890607533025901], [1.763328214900581, -2.894514914439559, -0.08123199086868012], [-1.011797100425901, -3.534969835652359, -0.743062865102722], [-0.821528296712382, -1.9466956411138796, -2.2051128650057947], [-2.011954187959717, -1.2102165034856942, 0.0009668077395988321], [-0.40994706586377744, -1.949605403411482, 1.259403781369355], [1.5770424436171686, 0.0, 0.0], [2.2927181468939204, 1.3915527243580526, 0.0], [2.341079856722333, 2.059852692894925, 1.416502376706476], [1.145391641007043, 2.033172555891592, 1.9759511228307955], [3.2022179922846545, 1.4422259961598203, 2.1984030628363276], [2.7215555591495058, 3.3152059779287892, 1.2745358845394041], [1.6292323391939847, 2.212255867310574, -0.809047933619886], [3.5455868300943827, 1.2600392214310698, -0.4280914688619856], [1.9974224573334847, -0.6906780683055271, 1.0535722235492986], [1.9277183224308958, -0.6529932317206277, -1.1102241252095293], [-0.3501493572534734, 0.5705349971623084, -1.1530217920585826], [-0.42668432219275476, 0.7576153073313041, 1.0049834283127295], [3.2753919637070847, -2.4456779116031013, -3.6069037652902414]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0016', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
Gaussian.command = 'g16 < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy()  # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    freq = reader_gauss.read_freq(logfile, ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'])
    zpe = reader_gauss.read_zpe(logfile)
    db.write(mol, name=label, data={'energy': e, 'frequencies': np.asarray(freq),
                                     'zpe': zpe, 'status': 'normal'})

except RuntimeError:
    for i in range(3):
        try:
            iowait(logfile, 'gauss')
            mol.positions = reader_gauss.read_geom(logfile, mol)
            kwargs = reader_gauss.correct_kwargs(logfile, kwargs)
            mol.calc = Gaussian(**kwargs)
            e = mol.get_potential_energy()  # use the Gaussian optimizer
            iowait(logfile, 'gauss')
            mol.positions = reader_gauss.read_geom(logfile, mol)
            freq = reader_gauss.read_freq(logfile, ['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'])
            zpe = reader_gauss.read_zpe(logfile)
            db.write(mol, name=label, data={'energy': e,
                                             'frequencies': np.asarray(freq),
                                             'zpe': zpe, 'status': 'normal'})
        except RuntimeError:
            if i == 2:
                db.write(mol, name=label, data={'status': 'error'})
            pass
        else:
            break

with open(logfile, 'a') as f:
    f.write('done\n')
