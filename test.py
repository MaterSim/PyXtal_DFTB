from pyxtal.db import database
from pyxtal.interface.dftb import DFTB
from pyxtal.elasticity import fit_elastic_constants, elastic_properties
import os
import warnings
warnings.filterwarnings("ignore")
import numpy as np
np.set_printoptions(precision=3)

#os.system('export OMP_NUM_THREADS=4')
skf_dir = os.environ['DFTB_PREFIX'] + '3ob-3-1/'
db_name, work_dir, codes = 'dataset/peptide.db', "peptide", ['IFABEW']
db_name, work_dir, codes = 'dataset/hydrocarbon.db', "PAH", ['BENZEN', 'NAPHTA']
db = database(db_name)

for code in codes:
    s = db.get_pyxtal(code)
    s.remove_water()
    symmetry = s.group.lattice_type
    for disp in ['TS', 'D3']:
        struc = s.to_ase(resort=True)
        for ftol in [1e-3]: #, 1e-4, 1e-5]:
            # Relax
            my = DFTB(struc, skf_dir, folder=work_dir, disp=disp)
            struc, energy = my.run(mode='vc-relax', step=3000, ftol=ftol)
            res = "\n{:8s} ".format(disp)
            for i in range(3): res += "{:7.4f} ".format(struc.cell[i, i])
            res += "{:12.4f}  ftol: {:.6f} ".format(energy, ftol)
            res += "Time: {:6.1f}".format(my.time)
            print(res)
    
            calc = my.get_calculator('relax', ftol=ftol)
            struc.set_calculator(calc)
            
            cwd = os.getcwd()
            os.chdir(work_dir)
            struc.write(code+'-'+disp+'.vasp', format='vasp', vasp5=True, direct=True)

            # Elastic
            C, C_err = fit_elastic_constants(struc, 
                                             symmetry=symmetry, 
                                             delta=1e-2, 
                                             verbose=False, 
                                             fmax=ftol)
            print(C)
            os.chdir(cwd)
