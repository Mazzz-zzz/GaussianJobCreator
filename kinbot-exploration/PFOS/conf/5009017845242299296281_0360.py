import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/home/akhalilov/GaussianJobCreator/kinbot-exploration/PFOS/kinbot.db')
label = 'conf/5009017845242299296281_0360'
logfile = 'conf/5009017845242299296281_0360.log'

mol = Atoms(symbols=['C', 'C', 'C', 'C', 'C', 'F', 'F', 'S', 'O', 'O', 'O', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'H'], positions=[[0.0, 0.0, 0.0], [-0.6937283445863813, 0.7718203945763873, 1.1635336229088487], [-2.2709622836291916, 0.7431123812655698, 1.1797556627389012], [-2.9699917885127123, -0.6604598778102058, 1.233811274632789], [-4.434530032795169, -0.5595652376518035, 1.7856565080301297], [-4.412149490303208, -0.463299714373019, 3.1035085495820747], [-5.049857983877613, 0.5028813412615539, 1.2743065103524127], [-5.456246996148601, -2.0667053949410255, 1.3314689650513938], [-4.653927843774284, -3.2240811782141816, 1.5294494737793916], [-6.738025277307821, -1.8824375218687037, 1.9045715392472673], [-5.582487214586102, -1.7982486697396678, -0.2153138328074454], [-3.0171454311974584, -1.1649147715918018, 0.004057906421360553], [-2.2839557356380094, -1.474510795144931, 2.025711167063385], [-2.62120547179292, 1.4189848017416544, 2.277641661887587], [-2.713377980741126, 1.3971249246107031, 0.11259346120335854], [-0.37104516182827885, 2.061737287215916, 1.0602591291106134], [-0.24552532002048744, 0.30567315029126324, 2.3240905646658536], [1.5770424436171646, 0.0, 0.0], [2.2927181468939195, 1.3915527243580557, 0.0], [2.341079856722324, 2.059852692894925, 1.4165023767064768], [1.1453916410070475, 2.0331725558915963, 1.975951122830796], [3.2022179922846545, 1.442225996159819, 2.1984030628363276], [2.721555559149513, 3.315205977928793, 1.2745358845394055], [1.6292323391939798, 2.2122558673105726, -0.8090479336198816], [3.545586830094381, 1.2600392214310707, -0.428091468861985], [1.9974224573334822, -0.6906780683055307, 1.053572223549296], [1.927718322430891, -0.6529932317206237, -1.1102241252095395], [-0.35014935725348145, -1.2838136616209426, 0.08241309473864378], [-0.4266843221927579, 0.4915335256355516, -1.1586058166012259], [-4.899930147116156, -2.285412412541228, -0.6977552531295398]])

kwargs = {'method': 'bmk', 'basis': '6-31G(d)', 'nprocshared': 8, 'mem': '700MW', 'label': 'conf/5009017845242299296281_0360', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'pop': 'None', 'freq': 'freq', 'opt': 'CalcFC'}
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
