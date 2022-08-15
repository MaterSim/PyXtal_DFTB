from pyxtal.db import database
from pyxtal.interface.dftb import DFTB
from pyxtal.elasticity import fit_elastic_constants
import ase.units as units
from ase.io import read
import os
import numpy as np
np.set_printoptions(formatter={'float': '{:8.3f}'.format})

skf_dir = '/home1/01606/qiangz/opt/dftb+/Dftb+sk/3ob-3-1/'
db_name, step, ncpu, ftol = 'dataset/peptide.db', 10000, 48, 1e-3
cmd = 'mpirun -np ' + str(ncpu) + ' /home1/01606/qiangz/opt/dftb+/bin/dftb+  > PREFIX.out'
os.environ['ASE_DFTB_COMMAND'] = cmd
db = database(db_name)

for code in ['IFABEW']:
    work_dir = code 
    s = db.get_pyxtal(code)
    s.remove_water()
    symmetry = s.group.lattice_type
    if symmetry == 'hexagonal': symmetry = 'trigonal_high'

    for disp in ['D3', 'TS']:
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
        C11, C33, C12, C13 = C[0, 0], C[2, 2], C[0, 1], C[0, 2]
        E1 = (C11**2*C33+2*C13**2*C12-2*C13**2*C11-C12**2*C33)/(C11*C33-C13**2)
        E3 = (C11**2*C33+2*C13**2*C12-2*C13**2*C11-C12**2*C33)/(C11**2-C12**2)

        strs = "\nYoung's Moduls (GPa): {:7.2f} {:7.2f}".format(E1, E3)
        print(strs)

        res += strs
        np.savetxt('Cij-'+disp+'.txt', C)
        with open('Status-'+disp+'.txt', 'w') as f:
            f.writelines(res)
        os.system("rm band.out detailed.out charges.bin results.tag geo_end.xyz")
        os.chdir(cwd)
