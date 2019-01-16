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

def processFile(file, max_waiting_time=max_waiting_time): 
    print file
    fig, ax = plt.subplots(1, 1, figsize=(60, 60))
    plt.xlim(-0.5, X+0.5)
    plt.ylim(-0.5, Y+0.5)
    ax.set_yticks(range(0, Y), minor=False)
    ax.set_xticks(range(0, X), minor=False)
    with open(file, 'r') as of:
        data = of.read().split('\n')
        xdata = []
        ydata = []
        cdata = []
        for line in data:
            if line != '':
                x, y, w = line.split(';')
                if float(w) >= 0:
                    xdata.append(int(x))
                    ydata.append(Y-int(y))
                    cdata.append(float(w))
    plt.scatter(xdata, ydata, c=cdata, edgecolors=None, s=(80.0), cmap='jet', vmin=-20.0, vmax=0.05*max_waiting_time)
    plt.savefig(file.split('_')[0]+"_py_out.jpg")
    plt.close(fig)

pool = multiprocessing.Pool(8)
pool.map(processFile, files)
#for file in files:

