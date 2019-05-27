from matplotlib import pyplot as plt
import pickle
from matplotlib.ticker import ScalarFormatter
from numpy import append

def plot(ct, ct2):
    fig, ax = plt.subplots()
    plt.loglog(ct[0][:], ct[1][:],'*-', markersize=15, label='hybrid')
    plt.loglog(ct2[0][:], ct2[1][:],'*-', markersize=15, label='python')
    plt.grid(True)
    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_major_formatter(ScalarFormatter())
        axis.set_minor_formatter(ScalarFormatter())
    yticks = ax.yaxis.get_minor_ticks()
    for t in yticks:
        t.label1.set_visible(False)
    xticks = ax.xaxis.get_minor_ticks()
    for t in xticks:
        t.label1.set_visible(False)
    plt.yticks(append(ct[1], ct2[1]))
    plt.xticks(ct[0])
    plt.legend(loc="upper left")
    plt.xlabel('Number of Grid Elements (log)',fontsize=20)
    plt.gcf().subplots_adjust(bottom=0.17, left=0.16)
    plt.ylabel('Runtime in s (log)',fontsize=20)
