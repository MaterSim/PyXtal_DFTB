from pyxtal.db import database
import os
import warnings; warnings.filterwarnings("ignore")
import numpy as np

db_name, ftol = '../dataset/hydrocarbon.db', 1e-3
db = database(db_name)

for code in db.codes:
    work_dir = code 
    if os.path.exists(work_dir):
        data = {}
        row = db.get_row(code)
        good=True
        for disp in ['TS', 'D3']:
            tag = "{:s}-{:s}".format(code, disp)
            F_cij = work_dir + '/Cij-' + tag + '.txt' 
            if os.path.exists(F_cij):
                C=np.loadtxt(F_cij)
                data['Cij_'+disp]  = C
            else:
                good=False
                break
        if good:
            db.db.update(row.id, data=data)
            print('write cij in ', code)
