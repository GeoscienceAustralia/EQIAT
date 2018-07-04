import os, sys
from glob import glob

for filename in glob(os.path.join('outputs', 'rupture*')):
    new_fn = filename.replace('.','p')
    new_fn = new_fn.replace('()','')
    new_fn = new_fn[:-4] + '.' + new_fn[-3:]
    os.rename(filename, new_fn)
#    cmd = 'mv ' + filename + ' ' + new_fn
#    print cmd
#    os.system(cmd)

