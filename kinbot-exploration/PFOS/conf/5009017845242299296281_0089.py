import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0089'
logfile = 'conf/5009017845242299296281_0089.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863827, 0.7718203945763839, 1.1635336229088482], [-2.270962283629194, 0.7431123812655638, 1.1797556627388983], [-3.020318448930585, 1.3845311280592592, -0.04033723253932561], [-4.4800581605680465, 0.8292457901518748, -0.18317744667967534], [-5.037241623014005, 0.7188730088199347, 1.010081276898237], [-5.21226269087264, 1.641135879597504, -0.9405485901848968], [-4.506227688342023, -0.8628257980321976, -0.994270676457445], [-5.768885931632611, -1.4622851118134734, -0.7320232749808816], [-3.9560519069161773, -0.7063724439402371, -2.2897187569198962], [-3.4184925635748487, -1.5627073800506124, -0.0956379932613879], [-3.086116569243988, 2.7007512384555885, 0.13918095834253777], [-2.363928874392461, 1.1178871860085964, -1.1621587584935842], [-2.596140664371263, -0.5518249635012428, 1.2269126406770372], [-2.689777026166652, 1.3325843308485112, 2.293046394029769], [-0.3710451618282821, 2.061737287215911, 1.0602591291106096], [-0.24552532002049118, 0.3056731502912633, 2.324090564665854], [1.5770424436171655, 0.0, 0.0], [2.2927181468939164, 1.3915527243580545, 0.0], [3.7823355744197134, 1.3186147352454611, -0.48079377647244326], [4.419592497958775, 0.349250732484395, 0.14976318527134813], [3.847332997716638, 1.1104432852891992, -1.7795740852228767], [4.3680104110515785, 2.46843408533521, -0.2056055442087179], [2.293896732220213, 1.8704189044736075, 1.2405689893126315], [1.6494649440008788, 2.2352148943143377, -0.8029305726284908], [1.9974224573334851, -0.6906780683055229, 1.0535722235492995], [1.927718322430896, -0.6529932317206226, -1.1102241252095322], [-0.350149357253474, -1.283813661620944, 0.08241309473865073], [-0.42668432219275637, 0.491533525635546, -1.1586058166012272], [-2.8240294132238324, -0.9027403877030197, 0.287793250356575]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0089', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
