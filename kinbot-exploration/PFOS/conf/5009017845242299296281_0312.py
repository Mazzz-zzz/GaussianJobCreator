import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0312'
logfile = 'conf/5009017845242299296281_0312.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863854, 0.6217394783082163, -1.2501828803165003], [-0.39761971585595896, -0.07566485901595352, -2.6339101198206674], [-1.1233882121466365, 0.4940238789133719, -3.902883790150537], [-1.2235971086413773, 2.0585875945462346, -3.8618506841213938], [-0.09347781129709946, 2.569072900004966, -3.4049731625056094], [-1.4669716767588572, 2.5390355666988684, -5.077925766590014], [-2.6155546136485244, 2.6391526295767536, -2.7451314376367266], [-2.591812817042398, 1.8537552012384486, -1.559796795240969], [-2.5907424079415757, 4.054970427660558, -2.759312013262729], [-3.829529432267254, 2.1746323601810986, -3.634529373821209], [-0.4343595853052552, 0.1444035564934193, -4.985471769970211], [-2.3533195674001606, 0.0018245094686737747, -3.976629593905572], [0.9186702629217385, 0.05926285465344133, -2.818308304723229], [-0.6863393184328452, -1.3666682794365117, -2.5232286857846664], [-0.2560445759534337, 1.8735740976390454, -1.3916701657561923], [-2.0076024771874494, 0.6427130616946748, -1.0543092166280594], [1.577042443617163, 0.0, 0.0], [2.292718146893913, 1.3915527243580548, 0.0], [3.7823355744197094, 1.3186147352454722, -0.4807937764724442], [4.4195924979587735, 0.34925073248439553, 0.1497631852713494], [3.847332997716632, 1.1104432852892097, -1.7795740852228779], [4.368010411051573, 2.4684340853352094, -0.20560554420871474], [2.2938967322202157, 1.870418904473611, 1.2405689893126353], [1.649464944000874, 2.23521489431434, -0.8029305726284839], [1.997422457333482, -0.6906780683055234, 1.0535722235493001], [1.9277183224308947, -0.6529932317206173, -1.110224125209532], [-0.3501493572534785, 0.7132786644586341, 1.0706086973199345], [-0.42668432219275504, -1.2491488329668508, 0.15362238828850236], [-4.5588438490409295, 2.8063796935101872, -3.564115064098955]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0312', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
