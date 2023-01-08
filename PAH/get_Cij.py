from pyxtal.db import database
import os
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from elastic import Elastic, ElasticOrtho, minimize, maximize

db_name, ftol = '../dataset/hydrocarbon.db', 1e-3
db = database(db_name)

count = 0
for code in db.codes:
    row = db.get_row(code)
    for key in ['Cij_D3', 'Cij_TS']:
        if key in row.data:
            success_load = True
            try:
                elas = Elastic(row.data[key])
                if elas.isOrthorhombic():
                    elas = ElasticOrtho(elas)
            except:
                success_load = False
                print(code, key)
                print("Cij", row.data[key])

            if success_load:
                minE = minimize(elas.Young, 2)[1]
                maxE = maximize(elas.Young, 2)[1]
                minLC = minimize(elas.LC, 2)[1]
                maxLC = maximize(elas.LC, 2)[1]
                #minG = minimize(elas.shear, 3)
			    #maxG = maximize(elas.shear, 3)
			    #minNu = minimize(elas.Poisson, 3)
			    #maxNu = maximize(elas.Poisson, 3)

                #print(minE)
                #print(maxE)
                print("{:12s} {:12s} {:12.3f} {:12.3f}".format(code, key, minLC, maxLC))
                count += 1
    #break

print(count)
