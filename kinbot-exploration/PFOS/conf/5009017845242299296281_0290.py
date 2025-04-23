import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0290'
logfile = 'conf/5009017845242299296281_0290.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586383, 0.6217394783082169, -1.2501828803165014], [-2.2709622836291934, 0.6501421835576519, -1.2334320314121738], [-3.0203184489305848, -0.7271986321270439, -1.1788705129599824], [-4.480058160568046, -0.5732592173009052, -0.6265591969129859], [-5.037241623014006, 0.5153195412709299, -1.1276029261821365], [-5.212262690872639, -1.6351069123925126, -0.950991067701113], [-4.5062276883420225, -0.4296507650339873, 1.2443643983651842], [-4.31336332463281, -1.7275122648777448, 1.792807533288349], [-3.719785844644731, 0.6988139336572083, 1.5811631839001288], [-6.018534089092507, -0.026581687886173513, 1.4204532938280494], [-3.086116569243988, -1.229841373580091, -2.408509660976093], [-2.363928874392461, -1.565402601090326, -0.3870393224017535], [-2.596140664371264, 1.3384499968011854, -0.1355618835040209], [-2.6897770261666523, 1.3195442638618278, -2.3005750802147817], [-0.3710451618282815, -0.11265730320380285, -2.3156464312138967], [-0.24552532002049296, 1.8598848945507218, -1.4267659957399788], [1.5770424436171662, 0.0, 0.0], [2.2927181468939173, 1.3915527243580565, 0.0], [1.6005215470082486, 2.4407219045638975, -0.9357086002340306], [1.3760692761371691, 1.9145350254105775, -2.1257143081021392], [0.46117091336275506, 2.8519467789832555, -0.4188289776134415], [2.404864907693503, 3.4780800111829135, -1.068930340330681], [3.540059223330468, 1.2295174124846107, -0.43152105569274457], [2.3088468039522536, 1.8960947387583758, 1.2310220414904738], [1.997422457333487, -0.6906780683055249, 1.0535722235492928], [1.9277183224308945, -0.6529932317206253, -1.1102241252095328], [-0.35014935725347374, 0.7132786644586341, 1.0706086973199322], [-0.42668432219275665, -1.2491488329668516, 0.15362238828850247], [-6.116585854580639, 0.5991948754171582, 2.151720639469122]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0290', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
