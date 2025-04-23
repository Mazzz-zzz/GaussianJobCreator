import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0440'
logfile = 'conf/5009017845242299296281_0440.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863817, 0.6217394783082193, -1.2501828803164992], [-2.2709622836291916, 0.6501421835576608, -1.2334320314121736], [-2.9970239643018983, 1.3180406141844534, -2.4534014845326175], [-3.0516123514917317, 2.8791655192442613, -2.3137795027818204], [-1.9012352660876388, 3.328390798231519, -1.8433210527854278], [-3.297380408862279, 3.4412767955173593, -3.493834169149421], [-4.411157127461294, 3.4273835647473945, -1.1421318134854248], [-4.152101609351178, 4.771999575916358, -0.758297176102659], [-5.635920736913452, 2.9825852620994855, -1.6966273883742395], [-4.063991433386536, 2.4809272846173105, 0.06792567294607363], [-2.3329497198238167, 1.017690771216737, -3.565960212041477], [-4.2414884028836, 0.8658956272092545, -2.5388589204157803], [-2.6427825570536188, -0.6321679360904127, -1.1859143708980966], [-2.6677893778920305, 1.2545909529572694, -0.12011753733594685], [-0.3710451618282831, -0.1126573032037984, -2.3156464312139007], [-0.24552532002048527, 1.859884894550721, -1.4267659957399743], [1.5770424436171664, 0.0, 0.0], [2.2927181468939204, 1.3915527243580523, 0.0], [3.7823355744197187, 1.3186147352454503, -0.48079377647244714], [4.419592497958777, 0.34925073248437954, 0.1497631852713498], [3.847332997716639, 1.1104432852891937, -1.779574085222872], [4.368010411051587, 2.468434085335197, -0.20560554420871657], [2.2938967322202246, 1.8704189044735993, 1.240568989312633], [1.6494649440008864, 2.235214894314337, -0.8029305726284844], [1.9974224573334791, -0.6906780683055327, 1.0535722235493006], [1.9277183224308905, -0.652993231720629, -1.1102241252095328], [-0.3501493572534734, 0.7132786644586355, 1.0706086973199367], [-0.4266843221927614, -1.2491488329668492, 0.15362238828849992], [-3.585329884639883, 1.6997836930134294, -0.24297655886264508]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0440', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
