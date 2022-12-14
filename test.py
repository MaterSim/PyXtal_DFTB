from pyxtal.db import database
from pyxtal.interface.dftb import DFTB
from pyxtal.elasticity import fit_elastic_constants, elastic_properties
from ase import units
from ase.io import read
import os
import numpy as np
import warnings; warnings.filterwarnings("ignore")

def run_calc(db, code, skf_dir, ftol, step1, step2, disps=['TS', 'D3']):
    """
    Workflow to relax structure and compute elastic constants

    Args:
        db:
        code:
        ftol:
        step1:
        step2:
        disps:
    """
    work_dir = code 
    s = db.get_pyxtal(code)
    s.remove_water()
    symmetry = s.group.lattice_type
    if symmetry == 'hexagonal': symmetry = 'trigonal_high'

    for disp in disps:
        tag = "{:s}-{:s}".format(code, disp)
        F_vasp = work_dir + '/Relaxed-' + tag + '.vasp'
        F_status = work_dir + '/Status-' + tag + '.txt' 
        F_cij = work_dir + '/Cij-' + tag + '.txt' 

        # Check if the results are converged
        converged = False
        if os.path.exists(F_status):
            with open(F_status, 'r') as f:
                line = f.readlines()[-1]
                if line.find('Converged') > -1:
                    converged = True

        if not converged:
            if os.path.exists(F_vasp):
                struc = read(F_vasp, format='vasp')
            else:
                struc = s.to_ase(resort=True)

            # Relax
            my = DFTB(struc, skf_dir, folder=work_dir, disp=disp)
            struc, energy = my.run(mode='vc-relax', step=step1, ftol=ftol)
            struc.write(F_vasp, format='vasp', vasp5=True, direct=True) 

            cwd = os.getcwd()
            os.chdir(work_dir)

            # Summary
            calc0 = my.get_calculator('single')
            struc.set_calculator(calc0)

            energy = struc.get_potential_energy()
            fmax = np.abs(struc.get_forces()).max()
            stress = struc.get_stress()/units.GPa

            res = "\n{:s} {:12.4f} eV {:6.1f} s\n".format(tag, energy, my.time)
            cell = struc.get_cell_lengths_and_angles()
            res += "Cell: {:7.4f} {:7.4f} {:7.4f} ".format(*cell[:3])
            res += "{:7.4f} {:7.4f} {:7.4f}\n".format(*cell[3:])
            res += "Stress (GPa) ".format(*stress)
            res += "{:7.2f} {:7.2f} {:7.2f} {:7.2f} {:7.2f} {:7.2f}\n".format(*stress)
            res += "Fmax  (eV/A) {:7.3f}\n".format(fmax)
            res += "ftol  (eV/A) {:.6f}\n".format(ftol)
            print(res)

            if fmax < ftol:
                finish = True
            else:
                finish = False

            # Elastic
            calc = my.get_calculator('relax', step=step2, ftol=ftol)
            struc.set_calculator(calc)

            C, C_err = fit_elastic_constants(struc, 
                                             symmetry=symmetry, 
                                             delta=1e-2, 
                                             verbose=False, 
                                             tag=tag,
                                             fmax=ftol)
            os.system("rm *.out *.hsd charges.bin geo_*")
            os.chdir(cwd)
            np.savetxt(F_cij, C)
        else:
            C = np.loadtxt(F_cij)
            res = 'Results from converged calculation\n'
            finish = True
            with open(F_status, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if line.find('Elastic Constants')>-1:
                        break
                    else:
                        res += line

        # Output Elastic constants
        res += "\nElastic Constants (GPa)\n"
        for i in range(6):
            strs = "{:7.2f} {:7.2f} {:7.2f} {:7.2f} {:7.2f} {:7.2f}".format(*C[i,:])
            res += strs + '\n'

        res += "\nBulk modulus, Shear modulus, Young's modulus, Poisson's ratio\n"
        k1, g1, e1, v1, k2, g2, e2, v2, k3, g3, e3, v3 = elastic_properties(C)
        res += "Voigt: {:7.2f} {:7.2f} {:7.2f} {:7.2f}\n".format(k1, g1, e1, v1)
        res += "Reuss: {:7.2f} {:7.2f} {:7.2f} {:7.2f}\n".format(k2, g2, e2, v2)
        res += "Hill:  {:7.2f} {:7.2f} {:7.2f} {:7.2f}\n".format(k3, g3, e3, v3)

  
        C11, C33, C12, C13 = C[0, 0], C[2, 2], C[0, 1], C[0, 2]
        E1 = (C11**2*C33+2*C13**2*C12-2*C13**2*C11-C12**2*C33)/(C11*C33-C13**2)
        E3 = (C11**2*C33+2*C13**2*C12-2*C13**2*C11-C12**2*C33)/(C11**2-C12**2)
        res += "\nYoung's Moduls (GPa): {:7.2f} {:7.2f}\n".format(E1, E3)

        print(res)
        with open(F_status, 'w') as f:
            f.write(res)
            if finish: f.write('Converged\n') 

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-n", "--ncpu", dest="ncpu",
                      help="number of cpus, default: 1",
                      type=int,
                      default=48,
                      metavar="ncpu")
    parser.add_option("-c", "--code", dest="code",
                      help="code, required",
                      metavar="code")
    parser.add_option("-d", "--db", dest="db",
                      help="database name, required",
                      metavar="code")
    parser.add_option("-s", "--s1", dest="s1",
                      type=int,
                      default=10000,
                      help="number of steps for Relaxation: 10000",
                      metavar="s1")
    parser.add_option("-e", "--s2", dest="s2",
                      type=int,
                      default=500,
                      help="number of steps for elastic: 500",
                      metavar="s2")
    parser.add_option("-f", "--ftol", dest="ftol",
                      type=float,
                      default=1e-3,
                      help="force tolerance, 1e-3 eV/A",
                      metavar="ftol")
    
    (options, args) = parser.parse_args()
    dbname = options.db
    code = options.code
    ncpu = options.ncpu
    ftol = options.ftol
    step1 = options.s1
    step2 = options.s2
    
    #Local
    #skf_dir = '/home/qzhu/opt/dftb+/Dftb+sk/3ob-3-1/'
    
    #Stampede2
    skf_dir = '/home1/01606/qiangz/opt/dftb+/Dftb+sk/3ob-3-1/'
    cmd = 'mpirun -np ' +str(ncpu)+ ' /home1/01606/qiangz/opt/dftb+/bin/dftb+ > PREFIX.out'
    os.environ['ASE_DFTB_COMMAND'] = cmd
    
    db = database(dbname)
    run_calc(db, code, skf_dir, ftol, step1, step2)
