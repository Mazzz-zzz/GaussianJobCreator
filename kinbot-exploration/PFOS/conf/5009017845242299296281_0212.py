import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0212'
logfile = 'conf/5009017845242299296281_0212.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863839, -1.3935598728845975, 0.08664925740765078], [-2.2709622836291916, -1.3932545648232184, 0.05367636867327128], [-2.9970239643018997, -2.783728318379927, 0.08524408716293168], [-3.051612351491735, -3.4433745877869124, -1.3365407299748326], [-1.9012352660876428, -3.2605582581586163, -1.9608104585981467], [-3.297380408862285, -4.746387544852171, -1.2333160417972255], [-4.411157127461302, -2.7028069473224603, -2.3971353288417934], [-5.653942072382576, -3.2473806773832288, -1.9712761490479032], [-4.145375387517529, -1.3154270127258754, -2.4966026875040415], [-4.014070723755674, -3.363990651540896, -3.7703436622992896], [-2.3329497198238167, -3.59705751812083, 0.9016340449500703], [-4.241488402883601, -2.631664135309425, 0.5195418500188208], [-2.642782557053618, -0.7109480038655852, 1.1404306775613267], [-2.6677893778920323, -0.7313203152515892, -1.0264488679511463], [-0.3710451618282831, -1.949079984012113, 1.255387302103286], [-0.24552532002048794, -2.16555804484198, -0.897324568925878], [1.5770424436171655, 0.0, 0.0], [2.292718146893915, 1.3915527243580577, 0.0], [1.6005215470082539, 2.4407219045638984, -0.9357086002340289], [1.3760692761371713, 1.914535025410575, -2.1257143081021477], [0.4611709133627504, 2.85194677898325, -0.41882897761345417], [2.404864907693501, 3.4780800111829064, -1.0689303403306945], [3.5400592233304655, 1.2295174124846056, -0.4315210556927474], [2.308846803952249, 1.8960947387583755, 1.2310220414904662], [1.997422457333483, -0.6906780683055255, 1.0535722235493015], [1.9277183224308954, -0.6529932317206226, -1.1102241252095353], [-0.35014935725347246, 0.5705349971623078, -1.1530217920585801], [-0.426684322192759, 0.7576153073313031, 1.0049834283127261], [-3.073932316052045, -3.592219815398517, -3.7740328114341604]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0212', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
