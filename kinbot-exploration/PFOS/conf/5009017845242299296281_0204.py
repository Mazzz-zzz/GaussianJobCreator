import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0204'
logfile = 'conf/5009017845242299296281_0204.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863848, 0.771820394576387, 1.1635336229088438], [-0.346602041513907, 0.2926212946843495, 2.625937724192378], [-0.9873929842445301, 1.0841115196807867, 3.8194659122849286], [-0.1847769512407539, 2.390826822806661, 4.147928663816875], [0.1883244129367716, 2.982413471416189, 3.026624712021102], [-0.9332368412513781, 3.2267902686432004, 4.861982548914763], [1.3577185656456106, 2.0276490900586026, 5.153151448047536], [2.218973759134343, 3.156231041734042, 5.069645293477194], [0.9276165596370575, 1.4463003806064307, 6.370691063047281], [1.9474063794356817, 0.8749010238476659, 4.256500563966419], [-2.23050689855434, 1.4246952716906223, 3.4913469454112596], [-1.0036314587762372, 0.32352185066805195, 4.906512125113516], [-0.7784534903451426, -0.9699349115591032, 2.6913872387857154], [0.9731582247379391, 0.3037494144073578, 2.768795082645787], [-2.0119541879597183, 0.6059455318059221, 1.0475948322279873], [-0.4099470658637772, 2.0654783699937878, 1.058705916025067], [1.5770424436171644, 0.0, 0.0], [2.2927181468939173, 1.3915527243580577, 0.0], [1.6005215470082557, 2.4407219045639006, -0.9357086002340307], [1.3760692761371713, 1.9145350254105749, -2.125714308102148], [0.46117091336275573, 2.851946778983254, -0.4188289776134514], [2.404864907693498, 3.4780800111829087, -1.0689303403306851], [3.5400592233304664, 1.2295174124846076, -0.4315210556927423], [2.308846803952246, 1.8960947387583749, 1.2310220414904718], [1.997422457333482, -0.6906780683055223, 1.053572223549301], [1.9277183224308967, -0.6529932317206268, -1.1102241252095277], [-0.3501493572534748, -1.2838136616209435, 0.0824130947386507], [-0.42668432219275526, 0.49153352563554525, -1.1586058166012319], [2.6003808887804003, 1.23101568359134, 3.637817263757279]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0204', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
