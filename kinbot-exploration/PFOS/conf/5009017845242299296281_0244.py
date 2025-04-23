import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0244'
logfile = 'conf/5009017845242299296281_0244.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863816, -1.393559872884599, 0.08664925740765088], [-0.346602041513903, -2.4204394252486647, -1.059551387211261], [-0.7363023803695204, -2.0273092602547935, -2.5275982746061083], [-0.5090926255515552, -0.4996509844589409, -2.8000337324856064], [-0.4446952607212927, -0.28040074331133114, -4.101680479585845], [-1.5014279159495583, 0.21624569992661208, -2.2788210606620134], [1.0867274347625948, 0.11649968183315958, -2.028425008687243], [1.4232745804972633, 1.3582849959413796, -2.6343427095121896], [0.9696038499690476, -0.07826862524410307, -0.6306753372183582], [2.052951777282571, -0.995549734914353, -2.585470109592334], [0.01837074119013288, -2.7227744393129876, -3.3736409548645963], [-2.0169579203546344, -2.300385488260373, -2.741562204025709], [0.9809271994195071, -2.5623271837262793, -1.0109564779824467], [-0.9145558447621683, -3.583976471922481, -0.7664486088494606], [-2.0119541879597156, -1.2102165034857004, 0.0009668077396009869], [-0.4099470658637734, -1.9496054034114851, 1.2594037813693522], [1.577042443617166, 0.0, 0.0], [2.292718146893914, 1.3915527243580579, 0.0], [2.3410798567223163, 2.059852692894934, 1.4165023767064617], [1.145391641007042, 2.0331725558915914, 1.9759511228307984], [3.202217992284646, 1.4422259961598303, 2.198403062836327], [2.7215555591494893, 3.3152059779287963, 1.2745358845394104], [1.629232339193968, 2.212255867310578, -0.8090479336198861], [3.5455868300943765, 1.2600392214310836, -0.4280914688619866], [1.997422457333485, -0.6906780683055282, 1.053572223549295], [1.927718322430894, -0.6529932317206228, -1.1102241252095346], [-0.35014935725347784, 0.5705349971623112, -1.1530217920585792], [-0.4266843221927563, 0.7576153073313007, 1.0049834283127304], [2.4864067372091454, -0.6896400595721657, -3.394480897495123]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0244', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
