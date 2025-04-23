import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0138'
logfile = 'conf/5009017845242299296281_0138.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.693728344586383, 0.6217394783082153, -1.250182880316502], [-0.397619715855956, -0.0756648590159588, -2.633910119820668], [1.0879414097563267, -0.08367235852192446, -3.138512306045848], [1.3412747433709995, -1.226837997824203, -4.181862535930469], [1.4942413294616161, -2.380495894559796, -3.555622467014096], [0.3200544491791535, -1.314882634219893, -5.029363151438022], [2.8849342358987173, -0.9110824602207875, -5.201213870219505], [2.5770387278690783, 0.07176458867049601, -6.1818478704616275], [3.9655135290297987, -0.8125637234609504, -4.291269244755091], [2.9826164228846865, -2.3155918755949196, -5.907155880784154], [1.350706773188926, 1.0832029916023946, -3.720189049227473], [1.9101736356829004, -0.2683335176786015, -2.1136883092367653], [-0.7843720924895166, -1.3445493546988787, -2.475221044816151], [-1.1581983763610533, 0.4930222734242909, -3.5614953002063583], [-0.25604457595342633, 1.8735740976390423, -1.3916701657561965], [-2.0076024771874454, 0.6427130616946737, -1.0543092166280619], [1.577042443617165, 0.0, 0.0], [2.2927181468939155, 1.3915527243580539, 0.0], [1.600521547008253, 2.4407219045638966, -0.9357086002340288], [1.3760692761371665, 1.9145350254105726, -2.1257143081021437], [0.4611709133627544, 2.851946778983249, -0.4188289776134503], [2.4048649076935025, 3.4780800111829073, -1.068930340330684], [3.5400592233304686, 1.2295174124846082, -0.43152105569274174], [2.3088468039522505, 1.8960947387583715, 1.2310220414904711], [1.9974224573334858, -0.6906780683055251, 1.0535722235493006], [1.927718322430896, -0.6529932317206267, -1.1102241252095306], [-0.3501493572534751, 0.7132786644586329, 1.0706086973199325], [-0.42668432219275654, -1.2491488329668492, 0.1536223882885044], [2.565677221337377, -2.2829517262051513, -6.779542669747005]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0138', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
