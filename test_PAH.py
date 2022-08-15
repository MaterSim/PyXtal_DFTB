from pyxtal.db import database
from pyxtal.interface.dftb import DFTB
from pyxtal.elasticity import fit_elastic_constants, elastic_properties
import ase.units as units
from ase.io import read
import os
import numpy as np
np.set_printoptions(formatter={'float': '{:8.3f}'.format})

skf_dir = '/home1/01606/qiangz/opt/dftb+/Dftb+sk/3ob-3-1/'
db_name, step, ftol = 'dataset/hydrocarbon.db', 1000, 1e-4
db = database(db_name)
 
cmd = 'mpirun -np 4 /home1/01606/qiangz/opt/dftb+/bin/dftb+  > PREFIX.out'
os.environ['ASE_DFTB_COMMAND'] = cmd
   
for code in ['BENZEN', 'NAPHTA']:
    work_dir = code
    s = db.get_pyxtal(code)
    s.remove_water()
    symmetry = s.group.lattice_type
    if symmetry == 'hexagonal': symmetry = 'trigonal_high'

    for disp in ['TS']: #, 'D3']:
        name = "{:s}-{:s}.vasp".format(code, disp)
        if os.path.exists(work_dir + '/' + name):
            struc = read(work_dir + '/' + name, format='vasp')
        else:
            struc = s.to_ase(resort=True)

        # Relax
        my = DFTB(struc, skf_dir, folder=work_dir, disp=disp)
        struc, energy = my.run(mode='vc-relax', step=step, ftol=ftol)

        cwd = os.getcwd()
        os.chdir(work_dir)
        
        # Summary
        calc0 = my.get_calculator('single')
        struc.set_calculator(calc0)

        res = "\n{:8s} ".format(disp)
        for i in range(3): res += "{:7.4f} ".format(struc.cell[i, i])
        res += "{:12.4f}  ftol: {:.6f} ".format(energy, ftol)
        res += "Time: {:6.1f}".format(my.time)
        stress = struc.get_stress()/units.GPa
        fmax = np.abs(struc.get_forces()).max()
        res += "\nStress (GPa) {:7.2f} {:7.2f} {:7.2f} {:7.2f} {:7.2f} {:7.2f}".format(*stress)
        res += "\nFmax  (eV/A) {:7.3f}".format(fmax)
        struc.write(name, format='vasp', vasp5=True, direct=True)
        print(res)

        # Elastic
        calc = my.get_calculator('relax', step=500, ftol=ftol)
        struc.set_calculator(calc)

        C, C_err = fit_elastic_constants(struc, 
                                         symmetry=symmetry, 
                                         delta=1e-2, 
                                         verbose=False, 
                                         fmax=ftol)
        print(C)
        np.savetxt('Cij-'+disp+'.txt', C)
        with open('Status-'+disp+'.txt', 'w') as f:
            f.writelines(res)
        os.system("rm band.out detailed.out charges.bin results.tag geo_end.xyz")
        os.chdir(cwd)
