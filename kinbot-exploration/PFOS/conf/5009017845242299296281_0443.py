import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0443'
logfile = 'conf/5009017845242299296281_0443.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863796, 0.6217394783082069, -1.250182880316509], [-2.2709622836291894, 0.6501421835576442, -1.2334320314121872], [-2.9699917885127105, 1.3987418462127608, -0.04493060495240165], [-4.434530032795166, 1.826206517213006, -0.40823054313394214], [-5.118760473256556, 2.0680039729419875, 0.6962442813079766], [-4.423212373540425, 2.9164461735590006, -1.169905961453703], [-5.334436835556408, 0.48700187924982774, -1.366542608756073], [-6.728895866558152, 0.7603575796655185, -1.3100476557255045], [-4.613374145520475, 0.28656904896842256, -2.568728970861194], [-5.018585672357864, -0.7262195873188438, -0.41319359314212145], [-3.017145431197457, 0.5859716358429903, 1.0068168322315652], [-2.283955735638007, 2.4915727289791856, 0.26410822321819105], [-2.6212054717929196, 1.2630031390416137, -2.367697716836102], [-2.713377980741125, -0.6010536646032328, -1.2662424075749643], [-0.3710451618282775, -0.11265730320381542, -2.315646431213901], [-0.2455253200204855, 1.8598848945507098, -1.426765995739992], [1.5770424436171655, 0.0, 0.0], [2.2927181468939155, 1.391552724358051, 0.0], [2.3410798567223225, 2.059852692894925, 1.4165023767064742], [1.1453916410070337, 2.0331725558916034, 1.9759511228307836], [3.2022179922846377, 1.4422259961598218, 2.198403062836331], [2.7215555591495093, 3.315205977928792, 1.2745358845394057], [1.6292323391939876, 2.2122558673105748, -0.8090479336198895], [3.545586830094388, 1.2600392214310632, -0.4280914688619763], [1.9974224573334811, -0.6906780683055206, 1.053572223549307], [1.9277183224308985, -0.6529932317206315, -1.1102241252095242], [-0.350149357253478, 0.7132786644586393, 1.0706086973199258], [-0.42668432219275826, -1.2491488329668514, 0.15362238828850694], [-4.214173792050014, -0.5531934929123324, 0.09567366279995665]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0443', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
