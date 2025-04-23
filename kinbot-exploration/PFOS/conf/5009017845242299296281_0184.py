import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0184'
logfile = 'conf/5009017845242299296281_0184.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863839, 0.6217394783082169, -1.2501828803164996], [-2.2709622836291916, 0.6501421835576581, -1.233432031412175], [-2.9970239643018997, 1.3180406141844505, -2.4534014845326175], [-3.051612351491735, 2.8791655192442587, -2.31377950278182], [-1.9012352660876428, 3.3283907982315184, -1.8433210527854291], [-3.297380408862285, 3.441276795517352, -3.4938341691494195], [-4.411157127461302, 3.427383564747387, -1.1421318134854237], [-5.653942072382576, 3.330865561641458, -1.826676087848644], [-4.145375387517529, 2.8198348568979403, 0.10910813390713721], [-4.014070723755674, 4.94720871831929, -1.0281295311781358], [-2.3329497198238167, 1.0176907712167336, -3.5659602120414777], [-4.241488402883601, 0.865895627209249, -2.5388589204157803], [-2.642782557053618, -0.6321679360904159, -1.1859143708980977], [-2.6677893778920323, 1.254590952957266, -0.1201175373359471], [-0.3710451618282831, -0.1126573032037984, -2.3156464312139007], [-0.24552532002048794, 1.8598848945507216, -1.4267659957399748], [1.5770424436171655, 0.0, 0.0], [2.2927181468939177, 1.3915527243580534, 0.0], [3.782335574419717, 1.3186147352454531, -0.48079377647244753], [4.419592497958776, 0.3492507324843861, 0.14976318527134874], [3.847332997716636, 1.1104432852891974, -1.779574085222872], [4.368010411051584, 2.4684340853352014, -0.20560554420871652], [2.293896732220222, 1.8704189044736026, 1.2405689893126326], [1.6494649440008828, 2.2352148943143386, -0.8029305726284839], [1.9974224573334807, -0.69067806830553, 1.0535722235492992], [1.9277183224308916, -0.6529932317206257, -1.1102241252095333], [-0.3501493572534748, 0.7132786644586343, 1.070608697319937], [-0.426684322192759, -1.2491488329668496, 0.15362238828849997], [-3.073932316052045, 5.0645181971172475, -1.2239372103958837]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0184', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
