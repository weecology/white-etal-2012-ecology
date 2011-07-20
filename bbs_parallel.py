import mete_sads
import numpy as np
import multiprocessing

pool = multiprocessing.Pool()

input_filename = '/home/kate/mete/bbs_dist_test.csv'

ifile = np.loadtxt(input_filename, delimiter = ',', 
                   dtype={'names': ('site','S','N'),
                          'formats':('S15', 'i8', 'f8')})

Niter = 1000
i=0
arg1=[];
while i < Niter:
    arg1.append(ifile)
    i+=1

pool.map(mete_sads.create_null_dataset, arg1)
pool.close()