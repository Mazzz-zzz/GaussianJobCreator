import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0441'
logfile = 'conf/5009017845242299296281_0441.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863835, 0.7718203945763832, 1.163533622908847], [-0.3466020415139052, 0.2926212946843425, 2.625937724192379], [-0.9873929842445267, 1.0841115196807767, 3.819465912284932], [-2.462172495902625, 1.515549635617316, 3.5056638996836136], [-3.0963703180309756, 1.7977976497721433, 4.630222090259054], [-2.4735048074690122, 2.582029157264965, 2.7110614524348216], [-3.414321394728243, 0.15733163429363328, 2.6279853329745646], [-4.8025540278360275, 0.44673458194602744, 2.7351632135037294], [-2.7471251593519987, -0.08691654408773557, 1.403016442551153], [-3.0690607739072786, -1.0296540093275808, 3.6038966073004945], [-0.9969507321948253, 0.3041218762439119, 4.896738560267487], [-0.2788854647228145, 2.1786256533403403, 4.065404334434962], [-0.7784534903451429, -0.9699349115591107, 2.691387238785714], [0.9731582247379417, 0.30374941440734815, 2.768795082645789], [-2.0119541879597196, 0.605945531805917, 1.047594832227991], [-0.4099470658637787, 2.065478369993785, 1.0587059160250714], [1.5770424436171657, 0.0, 0.0], [2.292718146893917, 1.3915527243580545, 0.0], [1.6005215470082557, 2.4407219045638984, -0.9357086002340258], [1.3760692761371698, 1.914535025410576, -2.125714308102144], [0.4611709133627575, 2.851946778983252, -0.4188289776134513], [2.4048649076935025, 3.4780800111829118, -1.0689303403306811], [3.5400592233304664, 1.2295174124846024, -0.43152105569274285], [2.3088468039522487, 1.8960947387583706, 1.2310220414904711], [1.997422457333483, -0.6906780683055259, 1.0535722235493021], [1.9277183224308978, -0.6529932317206234, -1.11022412520953], [-0.3501493572534738, -1.2838136616209443, 0.08241309473865076], [-0.4266843221927563, 0.49153352563554553, -1.158605816601229], [-3.7745605871978203, -1.144791638806447, 4.255798230569677]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0441', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
