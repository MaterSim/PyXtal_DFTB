from pyxtal.db import database
from pyxtal.interface.dftb import DFTB
from pyxtal.elasticity import fit_elastic_constants, elastic_properties
import ase.units as units
from ase.io import read
import os
import warnings; warnings.filterwarnings("ignore")
import numpy as np

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
    #print("+++++++++++++++++++++++++++++++++++++++++++++++++++", code)

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
            res += "Ftol: {:.6f}\n".format(ftol)
            cell = struc.get_cell_lengths_and_angles()
            res += "Cell: {:7.4f} {:7.4f} {:7.4f}  ".format(*cell[:3])
            res += "{:7.4f} {:7.4f} {:7.4f}\n".format(*cell[3:])
            res += "Stress (GPa) ".format(*stress)
            res += "{:7.2f} {:7.2f} {:7.2f} {:7.2f} {:7.2f} {:7.2f}\n".format(*stress)
            res += "Fmax  (eV/A) {:7.5f}\n".format(fmax)
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
            print(strs)

        res += "\nBulk modulus, Shear modulus, Young's modulus, Poisson's ratio\n"
        k1, g1, e1, v1, k2, g2, e2, v2, k3, g3, e3, v3 = elastic_properties(C)
        res += "Voigt: {:7.2f} {:7.2f} {:7.2f} {:7.2f}\n".format(k1, g1, e1, v1)
        res += "Reuss: {:7.2f} {:7.2f} {:7.2f} {:7.2f}\n".format(k2, g2, e2, v2)
        res += "Hill:  {:7.2f} {:7.2f} {:7.2f} {:7.2f}\n".format(k3, g3, e3, v3)

        print(res)
        with open(F_status, 'w') as f:
            f.write(res)
            if finish: f.write('Converged\n') 
#===============INPUTS================
if __name__ == "__main__":

    import multiprocessing as mp

    skf_dir = '/home1/01606/qiangz/opt/dftb+/Dftb+sk/3ob-3-1/'
    db_name, ftol = '../dataset/hydrocarbon.db', 1e-4
    db = database(db_name)
    
    nproc, total_cpu, step1, step2 = 48, 48, 5000, 500
    #nproc, total_cpu, step1, step2 = 48, 48, 200, 5#0#0
    ncpu = int(total_cpu/nproc)
    #cmd = 'mpirun -np '+str(ncpu)+ ' /home1/01606/qiangz/opt/dftb+/bin/dftb+ > PREFIX.out'
    cmd = '/home1/01606/qiangz/opt/dftb+/bin/dftb+ > PREFIX.out'
    os.environ['ASE_DFTB_COMMAND'] = cmd

    codes = [
             #'BCYBUE01', 'GESNIB', 'GOHWAB', 'TAVPEN', 'ZZZITY01','BENZEN', 'BOBLEJ', 
             #'HDMBZD', 'HAYXOU', 'HMBENZ04', 'KEGHEJ', 'NAPHTA', 'NPCBCP', 'SIWFUZ', 
             'XOZQOT', 'ANTCEN', 'BIPHEN', 'PHENAN', 'NPDCBU', 
             'CADWUD', 'NAFNIS', 'HIQHAT', 'TURTOR',
             'TAFZOP', 'SOYVIL', 'OROVAV', 'ZZZEMS01', 'NOFZIR', 'SUCSOY', 
             'DURENE', 'BUTPYR10', 'COPMUR', 'YAFJUL', 'MEYFUT', 
             'GUMZOE', 'HBZPEN', 'PYRCEN', 'DITBOX', 'KIDJUD', 
             'XITPEW', 'YETMEQ', 'CEKREP', 'CEKWEU', 'KIDJUD03', 
             'DIPRBZ', 'CURBIA', 'DIBENZ01', 'KOXQIA', 'DUPRIP', 
            ]
    #nproc = min([nproc, len(codes)])
    #pool = mp.Pool(processes=nproc)
    #for code in codes:
    #    pool.apply_async(run_calc, args=(db, code, skf_dir, ftol, step1, step2)) 
    #pool.close()
    #pool.join()

    queue = mp.Queue()
    per_cycle = int(np.ceil(total_cpu/ncpu))
    N_cycle = int(np.ceil(len(codes)/per_cycle))

    for cycle in range(N_cycle):
        N1, N2 = per_cycle * cycle, per_cycle * (cycle+1)
        if cycle + 1 == N_cycle: 
            N2 = len(codes)
    
        processes = []
        for i in range(N1, N2):
            p = mp.Process(target=run_calc,
                           args = (db, codes[i], skf_dir, ftol, step1, step2))
            p.start()
            processes.append(p)
        for p in processes: p.join()   
