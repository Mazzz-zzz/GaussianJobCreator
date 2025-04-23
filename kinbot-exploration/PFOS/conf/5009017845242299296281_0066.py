import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0066'
logfile = 'conf/5009017845242299296281_0066.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863835, -1.3935598728845968, 0.08664925740765074], [-0.3466020415139063, -2.420439425248663, -1.0595513872112614], [-0.9873929842445285, -3.8498102687678437, -0.9708648395635686], [-2.462172495902625, -3.793768812064689, -0.44032746470096895], [-3.0963703180309787, -4.908788780214294, -0.758172609562905], [-2.4735048074690127, -3.638862667661779, 0.8805721172861675], [-3.414321394728245, -2.354567876275693, -1.1777394743700764], [-3.1087457204169926, -2.283154722002885, -2.5648064073174446], [-4.729721703124768, -2.4221934601205444, -0.6576379643472848], [-2.667618736235436, -1.1834451528158607, -0.43514312627455315], [-0.9969507321948261, -4.392760927004429, -2.184992009459941], [-0.27888546472281445, -4.610056256946215, -0.14595700608828896], [-0.7784534903451441, -1.8458422644301227, -2.185681892820463], [0.9731582247379406, -2.549721586648359, -1.1213428320614884], [-2.01195418795972, -1.2102165034856947, 0.0009668077395977577], [-0.4099470658637775, -1.9496054034114818, 1.2594037813693502], [1.5770424436171657, 0.0, 0.0], [2.292718146893919, 1.3915527243580548, 0.0], [1.6005215470082612, 2.4407219045638997, -0.935708600234025], [1.376069276137171, 1.9145350254105789, -2.125714308102139], [0.46117091336276284, 2.8519467789832573, -0.4188289776134426], [2.404864907693509, 3.4780800111829064, -1.068930340330684], [3.5400592233304686, 1.2295174124846027, -0.43152105569274485], [2.3088468039522585, 1.8960947387583706, 1.2310220414904713], [1.9974224573334842, -0.6906780683055274, 1.053572223549297], [1.927718322430893, -0.6529932317206237, -1.1102241252095353], [-0.3501493572534762, 0.5705349971623119, -1.153021792058581], [-0.4266843221927546, 0.757615307331304, 1.0049834283127275], [-2.2748873847003046, -1.504313106331662, 0.38873078921492726]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0066', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
