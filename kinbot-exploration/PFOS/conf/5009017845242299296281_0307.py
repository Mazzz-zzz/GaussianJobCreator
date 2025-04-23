import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0307'
logfile = 'conf/5009017845242299296281_0307.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863787, 0.6217394783082173, -1.2501828803165034], [-2.27096228362919, 0.6501421835576563, -1.233432031412178], [-3.0203184489305843, -0.7271986321270375, -1.1788705129599861], [-2.2823447039971083, -1.7540796922695465, -0.2512417083665447], [-1.2544029753924617, -2.28103681980503, -0.8932338667375355], [-1.8529952461325283, -1.1533923865124285, 0.8549833339673658], [-3.4144078989664117, -3.1537286836453062, 0.27892708989652454], [-4.170179418850306, -3.5673580786500776, -0.852510553604157], [-2.639606570893646, -4.028839990459413, 1.0784393626921387], [-4.353747779728867, -2.334330219325272, 1.2416157077206607], [-4.24541884408069, -0.5314114745498586, -0.699508475208872], [-3.094058093800261, -1.247792535780131, -2.3970553321654844], [-2.596140664371261, 1.3384499968011907, -0.13556188350402626], [-2.689777026166647, 1.3195442638618329, -2.3005750802147857], [-0.371045161828279, -0.11265730320380046, -2.3156464312139007], [-0.24552532002048677, 1.8598848945507225, -1.4267659957399794], [1.5770424436171673, 0.0, 0.0], [2.292718146893919, 1.391552724358056, 0.0], [3.782335574419717, 1.318614735245452, -0.4807937764724453], [4.419592497958777, 0.3492507324843894, 0.1497631852713504], [3.8473329977166397, 1.110443285289194, -1.7795740852228756], [4.368010411051585, 2.4684340853351987, -0.20560554420871574], [2.2938967322202206, 1.870418904473604, 1.2405689893126335], [1.6494649440008837, 2.235214894314338, -0.8029305726284874], [1.997422457333483, -0.6906780683055304, 1.0535722235492997], [1.9277183224308945, -0.6529932317206277, -1.1102241252095342], [-0.3501493572534737, 0.713278664458636, 1.0706086973199311], [-0.42668432219275926, -1.2491488329668485, 0.15362238828849983], [-3.8941869095703203, -1.549722005458172, 1.5720035768190859]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0307', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
