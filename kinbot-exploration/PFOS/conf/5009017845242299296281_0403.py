import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0403'
logfile = 'conf/5009017845242299296281_0403.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863852, 0.621739478308219, -1.2501828803164985], [-0.3976197158559583, -0.07566485901595645, -2.6339101198206674], [-1.123388212146639, 0.4940238789133682, -3.9028837901505358], [-1.22359710864138, 2.058587594546232, -3.861850684121393], [-2.214848248727895, 2.4243645240243863, -3.0680431770923193], [-0.08476421785157562, 2.5840652181930435, -3.4193205106248867], [-1.5460987181777943, 2.781561421086144, -5.563310318993654], [-0.3251967235864591, 2.7648468452529986, -6.29238504106788], [-2.756206431121199, 2.212449625800765, -6.029258662900922], [-1.8531956580068718, 4.2632378624285385, -5.126424724669965], [-0.4343595853052585, 0.14440355649341485, -4.98547176997021], [-2.3533195674001606, 0.0018245094686675936, -3.976629593905572], [0.9186702629217383, 0.05926285465343887, -2.8183083047232307], [-0.6863393184328439, -1.3666682794365141, -2.5232286857846655], [-0.2560445759534327, 1.873574097639045, -1.3916701657561932], [-2.00760247718745, 0.6427130616946736, -1.0543092166280588], [1.5770424436171644, 0.0, 0.0], [2.292718146893911, 1.3915527243580579, 0.0], [3.7823355744197094, 1.3186147352454691, -0.4807937764724491], [4.419592497958774, 0.3492507324843974, 0.1497631852713452], [3.84733299771663, 1.110443285289208, -1.77957408522288], [4.368010411051577, 2.4684340853352156, -0.2056055442087264], [2.2938967322202144, 1.8704189044736128, 1.2405689893126337], [1.649464944000878, 2.235214894314343, -0.8029305726284891], [1.9974224573334887, -0.6906780683055254, 1.053572223549295], [1.927718322430894, -0.6529932317206204, -1.1102241252095357], [-0.3501493572534757, 0.7132786644586346, 1.0706086973199354], [-0.4266843221927576, -1.2491488329668534, 0.1536223882885027], [-2.1620326181532254, 4.2861951546864825, -4.20987953676282]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0403', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
