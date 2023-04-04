# Needs to have the access to CSD python API
# https://downloads.ccdc.cam.ac.uk/documentation/API/
from pyxtal.db import make_db_from_CSD

#codes = [
#         'IFABEW', #rigid
#         'CABJAU', #HYF
#         'ETIWIN', #0.75
#         'HIZCOJ', #2.5
#        ]
#codes = [
#         'BIYRIM01',
#         'DAHLOQ', 
#         'DAHMUX',
#         'YEWYAD',
#         'YEWYAD01',
#         'YEWYAD02',
#         'UHAWAD',
#         'UHAWIL',
#         'UHAWUX',
#         'UHAXEI',
#         'UCECAG01',
#         'UCECAG02',
#         'UCECAG03',
#         'VEWSIC',
#         'VEWSIC01',
#         'VEWSIC02',
#         'ACACCU06',
#         'DIBBUL01',
#         'HCLBNZ14',
#         'SOKDUU',
#         ]
#make_db_from_CSD('mech.db', codes)

codes = ['BENZEN',
         'UREAXX',
        ]
make_db_from_CSD('NPR.db', codes)
