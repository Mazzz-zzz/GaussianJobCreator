import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0469'
logfile = 'conf/5009017845242299296281_0469.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863845, 0.6217394783082141, -1.2501828803165025], [-0.397619715855956, -0.07566485901595821, -2.633910119820668], [-0.7491833517666574, -1.5988198371506703, -2.7680043658220965], [-2.081987459715323, -1.9555730430184117, -2.0223868238206664], [-2.976543962833394, -1.0033833036566044, -2.221397642773451], [-2.563171468208477, -3.1145119714526808, -2.463093699967502], [-1.8210733507287753, -2.1287237062650552, -0.17208792724638705], [-3.089985411108095, -2.050642212891716, 0.4652475191446673], [-0.9027569326874728, -3.1901545376136906, 0.01610431028305615], [-1.047509876275078, -0.7793973892984424, 0.07582470079310077], [-0.8953137364216059, -1.894503655492394, -4.056495166238402], [0.22590414679516443, -2.3330690032916617, -2.2478963639045055], [-1.13054897842193, 0.5952337900286009, -3.5269594944504683], [0.8878652548597804, 0.09037708117406053, -2.920716124729728], [-0.2560445759534272, 1.8735740976390407, -1.3916701657561954], [-2.0076024771874477, 0.6427130616946733, -1.0543092166280634], [1.5770424436171648, 0.0, 0.0], [2.292718146893916, 1.3915527243580563, 0.0], [3.7823355744197134, 1.3186147352454622, -0.48079377647244753], [4.4195924979587735, 0.34925073248439387, 0.14976318527135274], [3.8473329977166353, 1.1104432852891983, -1.7795740852228787], [4.3680104110515785, 2.4684340853352045, -0.20560554420871457], [2.293896732220217, 1.8704189044736097, 1.24056898931263], [1.6494649440008815, 2.2352148943143364, -0.8029305726284893], [1.9974224573334856, -0.6906780683055265, 1.0535722235492977], [1.9277183224308954, -0.6529932317206274, -1.1102241252095315], [-0.35014935725347507, 0.7132786644586356, 1.0706086973199327], [-0.42668432219275815, -1.2491488329668505, 0.15362238828850233], [-1.6603318471440074, -0.08897574638754503, 0.36518532996808084]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0469', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
