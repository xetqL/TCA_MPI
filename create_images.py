import os,sys
import matplotlib.pyplot as plt
import numpy as np

directory_content = os.listdir('.')
files = filter(lambda f: 'waiting_time' in f, directory_content)
files.sort()

X = int(sys.argv[1])
Y = int(sys.argv[2])
print X, Y
max_waiting_time = 0
max_wtime = []
x = 1/X
import multiprocessing
def getMaxFileSize(file):
    print "processing file: ", file
    global max_wtime
    max_waiting_time = 0
    with open(file, 'r') as of:
        data = of.read().split('\n')
        for line in data:
            if line != '':
                waiting_time = int(line.split(';')[-1])
                max_waiting_time = max(waiting_time, max_waiting_time)
    return max_waiting_time

pool = multiprocessing.Pool(8)
r = pool.map(getMaxFileSize, files)
pool.close()

max_waiting_time = max(r)
print max_waiting_time
import numpy as np
def processFile(file, max_waiting_time=max_waiting_time):
    print file
    arr = np.full((X+1, Y+1), 0)
    fig, ax = plt.subplots(1, 1, figsize=(20, 20))
    with open(file, 'r') as of:
        data = of.read().split('\n')
        for line in data:
            if line != '':
                x, y, w = line.split(';')
                if float(w) >= 0:
                    arr[int(x), Y-int(y)]= float(w)

    #ax[0].contourf(arr)
    ax.imshow(arr, cmap='jet', vmin=-5, vmax=15, interpolation='gaussian')

    m = plt.cm.ScalarMappable(cmap="jet")
    m.set_array(arr)
    m.set_clim(-5., 15)
    plt.savefig(file.split('_')[0]+"_py_out.jpg")
    plt.close(fig)

pool = multiprocessing.Pool(8)
pool.map(processFile, files)