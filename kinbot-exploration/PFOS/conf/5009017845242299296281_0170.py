import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0170'
logfile = 'conf/5009017845242299296281_0170.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863807, 0.7718203945763862, 1.1635336229088493], [-2.270962283629192, 0.7431123812655699, 1.1797556627389014], [-2.9699917885127123, -0.6604598778102058, 1.233811274632789], [-4.434530032795169, -0.559565237651802, 1.7856565080301303], [-5.030583914311105, 0.5082304477086564, 1.2846234874293212], [-5.127712885957739, -1.6471563109736436, 1.4607823325943794], [-4.4665010348623895, -0.41736338044220905, 3.656597399118879], [-4.226968110733983, -1.7075006466777913, 4.2048681537088655], [-3.7215317568045045, 0.7387617773290985, 3.993869049226321], [-5.992464116672816, -0.06933155880081121, 3.832248057772909], [-3.0171454311974584, -1.164914771591805, 0.0040579064213605645], [-2.2839557356380094, -1.4745107951449294, 2.0257111670633865], [-2.6212054717929205, 1.418984801741657, 2.2776416618875865], [-2.713377980741126, 1.397124924610703, 0.11259346120335914], [-0.37104516182827757, 2.061737287215914, 1.0602591291106138], [-0.24552532002048638, 0.3056731502912653, 2.324090564665853], [1.577042443617164, 0.0, 0.0], [2.292718146893918, 1.3915527243580552, 0.0], [2.341079856722327, 2.059852692894925, 1.4165023767064737], [1.145391641007045, 2.033172555891593, 1.9759511228307969], [3.2022179922846505, 1.4422259961598223, 2.198403062836329], [2.72155555914951, 3.315205977928792, 1.2745358845394055], [1.6292323391939774, 2.2122558673105734, -0.8090479336198829], [3.5455868300943783, 1.2600392214310698, -0.42809146886198174], [1.9974224573334816, -0.690678068305531, 1.0535722235492966], [1.9277183224308914, -0.6529932317206234, -1.11022412520954], [-0.35014935725347984, -1.2838136616209443, 0.08241309473864619], [-0.42668432219275787, 0.4915335256355502, -1.1586058166012263], [-6.113344565466942, 0.5523431323710007, 4.563596731219273]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0170', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
