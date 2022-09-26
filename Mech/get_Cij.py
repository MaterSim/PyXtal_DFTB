from pyxtal.db import database
import os
import warnings; warnings.filterwarnings("ignore")
import numpy as np

db_name = '../dataset/mech.db'
db = database(db_name)

count = 0
for code in db.codes:
    row = db.get_row(code)
    for key in ['Cij_D3', 'Cij_TS']:
        if key in row.data:
            print(code, row.data[key])
            count += 1
print(count)
