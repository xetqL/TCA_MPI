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

for file in files:
    with open(file, 'r') as of:
        data = of.read().split('\n')
        for line in data:
            if line != '':
                waiting_time = int(line.split(';')[-1])
                max_waiting_time = max(waiting_time, max_waiting_time)




print max_waiting_time


for file in files:
    A = np.eye(Y, X)
    A.fill(0.5)
    fig, ax = plt.subplots(1, 1, figsize=(3, 3))
    plt.xlim(-0.5, X+0.5)
    plt.ylim(-0.5, Y+0.5)
    ax.set_yticks(range(0, Y), minor=False)
    ax.set_xticks(range(0, X), minor=False)
    with open(file, 'r') as of:
        data = of.read().split('\n')
        for line in data:
            if line != '':
                x, y, w = line.split(';')
                plt.plot(int(x), Y-int(y),  's', color=(float(w) / max_waiting_time, 0, 0,  1))
    plt.grid(True)
    plt.savefig(file.split('_')[0]+"_py_out.jpg")
    plt.close(fig)