# Analyses and plots stuff

import matplotlib
matplotlib.use('TKAgg')


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os, sys

import re
import imp

def pause(message=''):
    print(message)
    programPause = raw_input("Press <ENTER> to continue...")

# To get parameters, taken from https://stackoverflow.com/questions/924700/best-way-to-retrieve-variable-values-from-a-text-file-python-json?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
def getVarFromFile(filename):
    import imp
    f = open(filename)
    global param
    param = imp.load_source('data', '', f)
    f.close()

# To sort filenames, taken from https://stackoverflow.com/questions/4623446/how-do-you-sort-files-numerically?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def getVector(filename):
    f = open(filename, 'r')
    a = f.readlines()
    a = [float(line.strip()) for line in a]
    f.close()
    return a

def plotMovie(a_files, b_files, c_files, s_files, filename):
    # Modified from:
    # Matplotlib Animation Example
    # author: Jake Vanderplas
    # email: vanderplas@astro.washington.edu
    # website: http://jakevdp.github.com
    # license: BSD
    fig = plt.figure()
    ax = plt.axes(xlim=(0, 0.6), ylim=(0, 10))
    line1, = ax.plot([], [], lw=2)
    line2, = ax.plot([], [], lw=2)
    line3, = ax.plot([], [], lw=2)
    line4, = ax.plot([], [], lw=2)

    # initialization function: plot the background of each frame
    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        line3.set_data([], [])
        line4.set_data([], [])
        return line1, line2, line3, line4

    # animation function.  This is called sequentially
    def animate(i):
        x = np.linspace(0,1,param.N_grid)
        a = getVector(path + a_files[i])
        b = getVector(path + b_files[i])
        c = getVector(path + c_files[i])
        s = getVector(path + s_files[i])
        # Renormalizing
        c = [c[j]*100 for j in range(0,len(c))]
        # s = [s[j] for j in range(0,len(c))]
        line1.set_data(x, a)
        line2.set_data(x, b)
        line3.set_data(x, c)
        line4.set_data(x, s)
        return line1, line2, line3, line4,

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=500, interval=300, blit=True)
    anim.save(filename, fps=15)
    plt.show()

'''
Stuff starts to happen here
'''

# Source
path='output/run1/'
# Get filenames
filenames = os.listdir(path)
# Sort to get correct time ordering
filenames = sorted(filenames, key=alphanum_key)
a_files = list()
b_files = list()
c_files = list()
s_files = list()
for file in filenames:
    if file[0]=="a":
        a_files.append(file)
    elif file[0]=="b":
        b_files.append(file)
    elif file[0]=="c":
        c_files.append(file)
    elif file[0]=="s":
        s_files.append(file)
    if file[0]=="p":
        p_file = file
# Get parameters
getVarFromFile(path + p_file)

# Plotting movie
plotMovie(a_files, b_files, c_files, s_files, 'diffusion_reaction_test.mp4')

# plt.plot(X,a)
# plt.show()
# pause('Showing plot.')
