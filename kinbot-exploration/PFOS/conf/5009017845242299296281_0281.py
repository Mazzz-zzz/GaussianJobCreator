import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0281'
logfile = 'conf/5009017845242299296281_0281.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863853, 0.7718203945763835, 1.1635336229088475], [-2.270962283629195, 0.7431123812655636, 1.1797556627388956], [-2.969991788512712, -0.6604598778102144, 1.233811274632784], [-3.0567398040787115, -1.3216748454235876, -0.185644503967384], [-3.2568593067482463, -2.6222799204493072, -0.06365857631325687], [-4.051038764929397, -0.7843632643516816, -0.8867620988473519], [-1.479977608212255, -1.0673537742871517, -1.170938927770445], [-0.3738206562783368, -1.2849959166013147, -0.30405963779567347], [-1.6585634249099073, -1.7284183817566876, -2.4104781415404664], [-1.6211807667533933, 0.47958742332824367, -1.4308951335960587], [-2.265942241823489, -1.4610539209006652, 2.0290187736421843], [-4.2028972120978745, -0.5310943056708419, 1.7067156695461185], [-2.6212054717929303, 1.4189848017416486, 2.2776416618875777], [-2.713377980741131, 1.3971249246106925, 0.11259346120335205], [-0.3710451618282853, 2.0617372872159114, 1.060259129110609], [-0.2455253200204914, 0.3056731502912629, 2.324090564665851], [1.5770424436171642, 0.0, 0.0], [2.2927181468939146, 1.3915527243580592, 0.0], [2.34107985672231, 2.059852692894934, 1.416502376706473], [1.1453916410070335, 2.0331725558915927, 1.9759511228307964], [3.202217992284636, 1.4422259961598374, 2.198403062836327], [2.721555559149494, 3.315205977928801, 1.2745358845394055], [1.6292323391939711, 2.212255867310578, -0.8090479336198867], [3.5455868300943783, 1.2600392214310814, -0.4280914688619822], [1.9974224573334831, -0.6906780683055265, 1.0535722235492986], [1.927718322430896, -0.6529932317206268, -1.1102241252095295], [-0.3501493572534738, -1.2838136616209457, 0.08241309473865084], [-0.426684322192755, 0.4915335256355444, -1.1586058166012285], [-1.1246355158056396, 0.9785806957138483, -0.767260509295052]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0281', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
