from pyxtal.db import database
from pyxtal.interface.dftb import DFTB
from ase.optimize import BFGS
from pyxtal.elasticity import fit_elastic_constants, elastic_properties
import os
import warnings
warnings.filterwarnings("ignore")
import numpy as np

np.set_printoptions(precision=3)
skf_dir = os.environ['DFTB_PREFIX'] + '3ob-3-1/'

work_dir = "tmp"
if not os.path.exists(work_dir): os.makedirs(work_dir)

db = database('dataset/hydrocarbon.db')

for code in ['BENZEN', 'NAPHTA']:
    s = db.get_pyxtal(code)
    symmetry = s.group.lattice_type
    print(s)
    for disp in ['TS', 'D3']:
        struc = s.to_ase()
        for ftol in [1e-3, 1e-4, 1e-5]:
            # Relax
            my = DFTB(struc, skf_dir, disp=disp)
            struc, energy = my.run(mode='vc-relax', step=2000, ftol=ftol)
            res = "\n{:8s} ".format(disp)
            res += "{:7.4f} ".format(struc.cell[0,0])
            res += "{:7.4f} ".format(struc.cell[1,1])
            res += "{:7.4f} ".format(struc.cell[2,2])
            res += "{:12.4f}  ftol: {:.6f} ".format(energy, ftol)
            res += "Time: {:6.1f}".format(my.time)
            print(res)
    
            calc = my.get_calculator('single')
            struc.set_calculator(calc)
    
            # Elastic
            C, C_err = fit_elastic_constants(struc, symmetry=symmetry, delta=1e-2, optimizer=BFGS, verbose=False, fmax=ftol)
            print(C)
