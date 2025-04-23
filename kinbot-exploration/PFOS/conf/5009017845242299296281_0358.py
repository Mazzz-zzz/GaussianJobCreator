import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0358'
logfile = 'conf/5009017845242299296281_0358.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863873, 0.6217394783082152, -1.2501828803164963], [-2.2709622836291974, 0.6501421835576506, -1.2334320314121674], [-2.9970239643019103, 1.3180406141844447, -2.453401484532607], [-4.47595427893113, 0.8159291580463225, -2.595649047209685], [-5.036363599701673, 0.7260147810600364, -1.402184032606544], [-5.17854859712389, 1.6536511332634092, -3.352935754365333], [-4.563626718048291, -0.8742349177524329, -3.4064031846728136], [-4.365065377906009, -0.7076207062110855, -4.80470361967359], [-3.808406466198901, -1.757900333220183, -2.597617319099752], [-6.087668062292554, -1.1732482806045788, -3.144808804612799], [-3.015076885801351, 2.635812688105754, -2.274118452082615], [-2.3510956819310325, 1.0276075788966117, -3.5753886090290967], [-2.642782557053619, -0.6321679360904244, -1.1859143708980922], [-2.6677893778920363, 1.2545909529572572, -0.12011753733593668], [-0.3710451618282891, -0.11265730320379726, -2.3156464312138985], [-0.24552532002049693, 1.8598848945507207, -1.4267659957399688], [1.5770424436171642, 0.0, 0.0], [2.2927181468939146, 1.3915527243580559, 0.0], [2.3410798567223283, 2.059852692894925, 1.4165023767064746], [1.1453916410070486, 2.033172555891597, 1.9759511228307982], [3.2022179922846528, 1.4422259961598274, 2.198403062836321], [2.721555559149504, 3.3152059779287955, 1.2745358845393997], [1.6292323391939727, 2.212255867310578, -0.8090479336198813], [3.5455868300943774, 1.260039221431077, -0.4280914688619914], [1.9974224573334856, -0.6906780683055329, 1.0535722235492941], [1.9277183224308876, -0.6529932317206254, -1.1102241252095377], [-0.35014935725347507, 0.7132786644586289, 1.070608697319939], [-0.42668432219275376, -1.249148832966853, 0.1536223882884981], [-6.6158401798865585, -0.9268757653497005, -3.9170116177514522]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0358', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
