import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0462'
logfile = 'conf/5009017845242299296281_0462.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863823, -1.3935598728845977, 0.08664925740765328], [-0.34660204151390084, -2.42043942524867, -1.0595513872112532], [-0.9873929842445218, -3.84981026876785, -0.9708648395635573], [-0.1847769512407465, -4.7876250073543885, -0.0034475673086342444], [0.889763588367667, -5.251777005738083, -0.6169844400485547], [0.18128897051813278, -4.122301251114864, 1.0884499008005055], [-1.211921939867945, -6.25481490687108, 0.5569008084445429], [-2.1177150707628143, -5.807327999797046, 1.5577630766068387], [-1.6023819083198627, -6.953526781881375, -0.6113027270766194], [-0.06536911483788746, -7.084547967303732, 1.2479162190177586], [-2.2305068985543333, -3.735942783996665, -0.5118511747699884], [-1.0036314587762272, -4.410925069658706, -2.173077921198876], [-0.7784534903451336, -1.845842264430133, -2.1856818928204595], [0.9731582247379497, -2.5497215866483645, -1.121342832061476], [-2.0119541879597156, -1.210216503485697, 0.0009668077395988345], [-0.4099470658637802, -1.94960540341148, 1.259403781369354], [1.5770424436171675, 0.0, 0.0], [2.2927181468939195, 1.3915527243580532, 0.0], [3.782335574419716, 1.3186147352454565, -0.48079377647244614], [4.419592497958777, 0.34925073248438965, 0.14976318527135635], [3.847332997716641, 1.11044328528919, -1.7795740852228872], [4.368010411051579, 2.4684340853351987, -0.20560554420871924], [2.2938967322202135, 1.870418904473607, 1.2405689893126295], [1.6494649440008882, 2.235214894314335, -0.8029305726284909], [1.9974224573334811, -0.6906780683055251, 1.0535722235493044], [1.9277183224308978, -0.6529932317206305, -1.1102241252095304], [-0.3501493572534712, 0.5705349971623052, -1.153021792058584], [-0.4266843221927554, 0.7576153073313081, 1.0049834283127217], [-0.04872451318222775, -6.9074675804434715, 2.1988777747712738]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0462', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
