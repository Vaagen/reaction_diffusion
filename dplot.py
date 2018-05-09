# Analyses and plots stuff

import matplotlib
matplotlib.use('TKAgg')


import numpy as np
import peakutils
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os, sys

import re
import imp

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)

def pause(message=''):
    print(message)
    programPause = input("Press <ENTER> to continue...")

def getVarFromFile(filename):
    # To get parameters, taken from https://stackoverflow.com/questions/924700/best-way-to-retrieve-variable-values-from-a-text-file-python-json?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
    import imp
    f = open(filename)
    global param
    param = imp.load_source('parameters', '', f)
    f.close()

def tryint(s):
    # To sort filenames, taken from https://stackoverflow.com/questions/4623446/how-do-you-sort-files-numerically?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    # To sort filenames, taken from https://stackoverflow.com/questions/4623446/how-do-you-sort-files-numerically?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
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

def plotMovie(a_files, b_files, c_files, s_files, filename, legend=None):
    # Modified from:
    # Matplotlib Animation Example
    # author: Jake Vanderplas
    # email: vanderplas@astro.washington.edu
    # website: http://jakevdp.github.com
    # license: BSD
    fig = plt.figure()
    ax = plt.axes(xlim=(0, 1), ylim=(0, 10))
    line1, = ax.plot([], [], lw=2)
    line2, = ax.plot([], [], lw=2)
    line3, = ax.plot([], [], lw=2)
    line4, = ax.plot([], [], lw=2)
    addlabel()
    fig.legend(legend,fontsize=14,loc=0)

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
        addlabel()
        fig.legend(legend,fontsize=14,loc=0)
        return line1, line2, line3, line4,

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=500, interval=300, blit=True)
    # anim.save(filename, fps=15)
    plt.show()

def addParameters(fig, ptxt):
    fig.text(1, .1, ptxt, position=(0.85,0.2),fontsize=14)

def addlabel():
    plt.xlabel('X/L')
    plt.ylabel('Density')


if __name__ == "__main__":
    # Source
    path='output/run2/'
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
        if file[0]=="t":
            t_nucleation_file= file
    # Get parameters
    getVarFromFile(path + p_file)
    # Make text for parameters
    param_names = [item for item in dir(param) if not item.startswith("__") and not item.startswith("CFL")]
    param_vals = [getattr(param, name) for name in param_names]
    param_text = [ param_names[i] + '= '+ '{0:.2e}'.format(param_vals[i]) for i in range(0,len(param_vals)) if not param_vals[i] == None]
    ptxt=''
    for i in range(0,len(param_text)):
        ptxt = ptxt + param_text[i] + '\n'
    # For plots with all functions
    legend = ['a', 'b', '100*c', 's']

    # Plotting movie
    # plotMovie(a_files, b_files, c_files, s_files, 'diffusion_reaction_test.mp4',legend=legend)

    # Looking at peaks
    # At frame frame i
    frame = 200
    x = np.linspace(0,1,param.N_grid)
    s = getVector(path + s_files[frame])
    # Finding peaks
    peaks = peakutils.indexes(s, thres=0.1/max(s), min_dist=0.01)
    peaknums = range(1,len(peaks))
    peak_x = [x[j] for j in peaks]
    peak_s = [s[j] for j in peaks]
    fig1=plt.figure(1)
    plt.plot(x,s)
    plt.plot(peak_x,peak_s,'rx')
    # Space between peaks
    fig2=plt.figure(2)
    delta_peak_x = np.diff(peak_x)
    plt.plot(peaknums, list(reversed(delta_peak_x)))
    # Time between nucleation
    fig3 = plt.figure(3)
    peak_t = getVector(path + t_nucleation_file)
    peak_t = [peak_t[j] for j in peaks]
    delta_peak_t = -np.diff(peak_t)
    plt.plot(peaknums, list(reversed(delta_peak_t*param.delta_t)))

    addParameters(fig1,ptxt)
    addlabel()
    plt.legend(legend)

    plt.show()
    # plt.show(block=False)
    # pause('Showing plot.')
