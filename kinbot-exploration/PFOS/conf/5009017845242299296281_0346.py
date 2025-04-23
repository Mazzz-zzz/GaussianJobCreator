import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0346'
logfile = 'conf/5009017845242299296281_0346.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863841, -1.393559872884598, 0.0866492574076533], [-0.34660204151390467, -2.4204394252486665, -1.0595513872112565], [1.1624535490467052, -2.799902192493112, -1.2603960465750403], [1.3229555267971602, -4.171696370783135, -2.0033664407024236], [0.435797104589275, -4.260501353184641, -2.978848655857206], [2.548490142211768, -4.284594543793278, -2.5075677392970896], [1.0607552256216175, -5.626043159278572, -0.8467378261629684], [0.824409821566643, -6.785122611371573, -1.6361088932333756], [2.05686193726949, -5.53870203229828, 0.15600826243368748], [-0.30216120032047983, -5.164592233814887, -0.2063750468769337], [1.750051457246917, -1.8550381832545373, -1.989060753302588], [1.7633282149005816, -2.8945149144395605, -0.08123199086867887], [-1.011797100425899, -3.534969835652361, -0.7430628651027225], [-0.8215282967123796, -1.9466956411138807, -2.205112865005796], [-2.011954187959719, -1.2102165034856955, 0.0009668077395988332], [-0.40994706586378016, -1.9496054034114823, 1.2594037813693555], [1.5770424436171662, 0.0, 0.0], [2.2927181468939177, 1.3915527243580543, 0.0], [3.782335574419718, 1.3186147352454556, -0.4807937764724424], [4.419592497958777, 0.349250732484388, 0.14976318527135196], [3.8473329977166397, 1.1104432852891937, -1.7795740852228779], [4.368010411051582, 2.4684340853351987, -0.20560554420871624], [2.293896732220218, 1.8704189044736048, 1.2405689893126308], [1.6494649440008855, 2.235214894314336, -0.8029305726284819], [1.9974224573334822, -0.6906780683055258, 1.0535722235492992], [1.927718322430895, -0.6529932317206273, -1.1102241252095286], [-0.35014935725347296, 0.5705349971623079, -1.1530217920585832], [-0.4266843221927559, 0.7576153073313044, 1.0049834283127261], [-1.04656261460974, -5.574276546977027, -0.6689630759579035]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0346', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
