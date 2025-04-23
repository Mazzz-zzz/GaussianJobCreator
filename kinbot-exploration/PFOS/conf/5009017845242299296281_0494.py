import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0494'
logfile = 'conf/5009017845242299296281_0494.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863837, 0.7718203945763864, 1.1635336229088453], [-2.2709622836291934, 0.7431123812655719, 1.1797556627388928], [-2.969991788512716, -0.6604598778102042, 1.233811274632781], [-3.0567398040787155, -1.3216748454235774, -0.18564450396738677], [-4.057927071959131, -0.792630756448847, -0.866994890796401], [-1.9208515428855029, -1.1373859158669193, -0.8524994275297271], [-3.34387825427727, -3.1729541587223338, -0.0761924793873893], [-2.111023082257976, -3.7949253398575564, 0.26388248657236163], [-4.551026485388037, -3.3601655975089146, 0.6401187636917787], [-3.6449190406827388, -3.4468180380652536, -1.597677165580655], [-2.2659422418234954, -1.461053920900658, 2.0290187736421816], [-4.202897212097876, -0.531094305670829, 1.7067156695461143], [-2.6212054717929294, 1.4189848017416595, 2.2776416618875723], [-2.713377980741128, 1.3971249246107003, 0.11259346120334827], [-0.37104516182828085, 2.061737287215915, 1.060259129110606], [-0.2455253200204934, 0.30567315029126696, 2.3240905646658496], [1.577042443617164, 0.0, 0.0], [2.2927181468939195, 1.3915527243580508, 0.0], [1.600521547008264, 2.4407219045638957, -0.9357086002340329], [1.376069276137176, 1.9145350254105789, -2.125714308102145], [0.4611709133627717, 2.8519467789832573, -0.4188289776134594], [2.404864907693521, 3.4780800111829038, -1.0689303403306856], [3.54005922333047, 1.2295174124845951, -0.43152105569274307], [2.308846803952257, 1.8960947387583684, 1.231022041490474], [1.9974224573334791, -0.6906780683055299, 1.0535722235492997], [1.927718322430897, -0.6529932317206304, -1.1102241252095253], [-0.35014935725347784, -1.2838136616209423, 0.08241309473865063], [-0.42668432219275393, 0.4915335256355461, -1.15860581660123], [-2.8416156026699833, -3.736555956700702, -2.052339580651094]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0494', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
