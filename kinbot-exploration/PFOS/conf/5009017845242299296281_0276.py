import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0276'
logfile = 'conf/5009017845242299296281_0276.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863829, -1.3935598728845986, 0.08664925740765582], [-2.2709622836291947, -1.3932545648232149, 0.05367636867327487], [-2.997023964301903, -2.7837283183799224, 0.08524408716294145], [-3.051612351491736, -3.443374587786916, -1.33654072997482], [-3.2793160186248658, -4.740067389005456, -1.2214722185603628], [-4.013520036470854, -2.8879361664696814, -2.0680609321946553], [-1.4404858172046033, -3.2181908870364606, -2.2720763640087394], [-1.3958302980578645, -1.8855634219156592, -2.7665871643209754], [-0.4216488399427635, -3.806435541107727, -1.4838529106856715], [-1.736726766178271, -4.192925703027661, -3.4731908333474366], [-2.332949719823824, -3.597057518120824, 0.9016340449500823], [-4.241488402883606, -2.6316641353094186, 0.5195418500188294], [-2.6427825570536196, -0.7109480038655773, 1.1404306775613278], [-2.6677893778920323, -0.7313203152515901, -1.0264488679511437], [-0.3710451618282885, -1.949079984012104, 1.255387302103295], [-0.24552532002048857, -2.165558044841984, -0.8973245689258686], [1.5770424436171655, 0.0, 0.0], [2.2927181468939173, 1.3915527243580526, 0.0], [1.6005215470082599, 2.4407219045638957, -0.9357086002340331], [1.3760692761371762, 1.9145350254105757, -2.1257143081021517], [0.4611709133627524, 2.8519467789832524, -0.4188289776134535], [2.404864907693506, 3.478080011182909, -1.0689303403306822], [3.5400592233304717, 1.229517412484606, -0.43152105569273425], [2.308846803952246, 1.8960947387583724, 1.2310220414904702], [1.9974224573334807, -0.6906780683055245, 1.0535722235493075], [1.9277183224308958, -0.6529932317206303, -1.1102241252095275], [-0.35014935725347246, 0.5705349971623062, -1.1530217920585848], [-0.4266843221927599, 0.757615307331308, 1.0049834283127232], [-0.9262218000357515, -4.644489097915783, -3.7472977229934794]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0276', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
