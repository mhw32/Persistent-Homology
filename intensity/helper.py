import numpy as np
# import sys
# sys.path.append('../correlation')
from os import listdir
from os.path import isfile, join
# import translate

def apply_to_foam(foam, func):
    return [apply_to_vector(vec, func) for vec in foam]

def apply_to_vector(vector, func):
    return np.array([func(i) for i in vector])

def files_in_folder(folder):
    onlyfiles = [f for f in listdir(folder) if isfile(join(folder, f))]
    return onlyfiles

# def convert_folder_rds_to_np(rdsfolder, npfolder):
#     rdsfiles = files_in_folder(rdsfolder)
#     filecount = 0
#     for files in rdsfiles:
#         if files[-4:] == '.rds':
#             filecount += 1
#             print('converting %s' % files)
#             tmp = translate.rds_to_np(join(rdsfolder, files))
#             simplify = lambda x: [np.array(i) for i in x]
#             if type(tmp[0][0]) == np.float:
#                 tmp = np.array(simplify(tmp))
#             else:
#                 tmp = np.array([simplify(i) for i in tmp])
#             newfiles = files[:-4] + '.npy'
#             np.save(join(npfolder, newfiles), tmp)
    
#     print("finished. %d files converted." % filecount)

