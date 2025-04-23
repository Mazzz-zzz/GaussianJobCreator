import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0457'
logfile = 'conf/5009017845242299296281_0457.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863828, -1.3935598728845988, 0.08664925740764715], [-2.2709622836291925, -1.393254564823215, 0.05367636867326867], [-3.0203184489305857, -0.6573324959322177, 1.2192077454993016], [-2.2823447039971074, 0.6594581441991476, 1.6446984279510903], [-1.8454205729265487, 1.2959521938724365, 0.5720966511971384], [-3.1051738563474744, 1.4520556038051013, 2.325610678994924], [-0.8077838749090341, 0.3063105249128043, 2.750420116002146], [0.02404031974969211, 1.4596991880381198, 2.764352408145413], [-1.308885263300432, -0.32452335462936927, 3.9150242490113], [-0.12061752729807175, -0.8013021432139203, 1.866373048292962], [-4.2454188440806915, -0.34008637241846285, 0.8099700744271566], [-3.0940580938002635, -1.452014544042186, 2.279147700720935], [-2.596140664371264, -0.7866250332999404, -1.091350757173021], [-2.689777026166652, -2.652128594710337, 0.007528686185003316], [-0.37104516182828423, -1.9490799840121118, 1.2553873021032829], [-0.24552532002048819, -2.1655580448419807, -0.8973245689258816], [1.5770424436171657, 0.0, 0.0], [2.2927181468939155, 1.3915527243580548, 0.0], [1.6005215470082557, 2.4407219045638957, -0.9357086002340326], [1.3760692761371658, 1.9145350254105729, -2.1257143081021472], [0.4611709133627617, 2.8519467789832533, -0.4188289776134545], [2.4048649076935122, 3.478080011182908, -1.0689303403306878], [3.540059223330465, 1.229517412484603, -0.43152105569275057], [2.3088468039522545, 1.8960947387583758, 1.2310220414904705], [1.9974224573334807, -0.6906780683055267, 1.0535722235493008], [1.9277183224308954, -0.6529932317206231, -1.110224125209533], [-0.3501493572534757, 0.5705349971623115, -1.1530217920585826], [-0.42668432219275804, 0.7576153073313027, 1.004983428312724], [0.5653126388898282, -0.40842171463951105, 1.3086032367643816]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0457', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
