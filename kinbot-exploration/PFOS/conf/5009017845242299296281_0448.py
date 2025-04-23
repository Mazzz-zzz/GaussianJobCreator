import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0448'
logfile = 'conf/5009017845242299296281_0448.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863833, 0.771820394576383, 1.1635336229088467], [-0.39761971585595657, 2.318865504557593, 1.2514273698287435], [1.0879414097563251, 2.759867566386748, 1.4967937649483776], [1.3412747433709984, 4.235018190162304, 1.028458395521435], [1.4942413294616166, 4.269507329980803, -0.2837586847863007], [0.3200544491791523, 5.0129975711126376, 1.3759598114895801], [2.884934235898715, 4.959924572236465, 1.8115863796161271], [3.2682214306983735, 6.105637772983566, 1.0613442881906772], [2.6730330558510453, 4.962960223553066, 3.21174545856519], [3.8803541094764227, 3.7845843027522963, 1.4825618273036472], [1.3507067731889246, 2.6801767277104727, 2.7981758327967103], [1.9101736356829004, 1.964674530320517, 0.8244605116218737], [-0.7843720924895166, 2.8158789821420864, 0.07319662459686944], [-1.1581983763610555, 2.8378342687254476, 2.207717463520169], [-0.2560445759534273, 0.2684346684142449, 2.3183978473060134], [-2.0076024771874477, 0.5917020341966364, 1.0837604470856907], [1.5770424436171657, 0.0, 0.0], [2.2927181468939146, 1.3915527243580563, 0.0], [1.6005215470082461, 2.4407219045638993, -0.9357086002340281], [1.3760692761371636, 1.9145350254105729, -2.125714308102147], [0.4611709133627515, 2.851946778983249, -0.4188289776134535], [2.40486490769349, 3.4780800111829135, -1.068930340330684], [3.540059223330462, 1.2295174124846096, -0.4315210556927445], [2.3088468039522456, 1.8960947387583755, 1.2310220414904718], [1.997422457333485, -0.6906780683055237, 1.0535722235493017], [1.9277183224308967, -0.6529932317206261, -1.1102241252095297], [-0.3501493572534727, -1.2838136616209468, 0.08241309473865208], [-0.4266843221927559, 0.491533525635541, -1.1586058166012299], [3.39476229903792, 2.9563490114953845, 1.363443051925309]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0448', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
